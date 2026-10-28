#!/usr/bin/env python3
"""PARMELA orbit correction: sequential one-to-one or coupled SVD.

Usage: python autocorrection.py input.inp -1vs1|-svd delta_val iterations
Standard library only. See autocorrection_README.md.
"""

import argparse
from bisect import bisect_left
import csv
from dataclasses import dataclass
from datetime import datetime
import hashlib
import math
from pathlib import Path
import re
import subprocess
import sys
import time

# Acceptance is per plane: abs(<X>) < tolerance AND abs(<Y>) < tolerance.
BPM_TOLERANCE_MM = 0.03
# Safeguard the predicted vertex; enlarge the scan on later sweeps if needed.
MAX_STEP_MULTIPLIER = 4.0
# Discard singular modes weaker than this fraction of the largest response.
SVD_RCOND = 1e-3
SVD_BACKTRACK_STEPS = 6
# Relative error in the predicted orbit CHANGE before a matrix is rebuilt.
SVD_PREDICTION_TOLERANCE = 0.25

# Other named elements are discovered from OUTPAR. POISSON and CFIELD describe
# fields and do not occupy element indices in this PARMELA version.
ELEMENT_NAMES = {"drift", "cathode", "solenoid", "steerer", "quad",
                 "cell", "trwave", "bend"}
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
ELEMENT_HEADER = re.compile(r"^\s*n\s+z1\s+element\s+z2\b", re.I)
ELEMENT_ROW = re.compile(
    r"^\s*(\d+)\s+(" + NUMBER + r")\s+(.*?)\s+(" + NUMBER
    + r")\s+(" + NUMBER + r")(?:\s|$)")


def timestamp():
    """Local wall-clock time with an explicit UTC offset."""
    return datetime.now().astimezone().isoformat(timespec="seconds")


def elapsed_text(seconds):
    """Readable elapsed time; hours may exceed 24 for long corrections."""
    milliseconds = round(seconds * 1000)
    hours, remainder = divmod(milliseconds, 3600000)
    minutes, remainder = divmod(remainder, 60000)
    whole_seconds, fraction = divmod(remainder, 1000)
    return "%02d:%02d:%02d.%03d" % (hours, minutes, whole_seconds, fraction)


def number(value):
    result = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError("Non-finite number: %s" % value)
    return result


def code_part(line):
    return re.split(r"[!;]", line, maxsplit=1)[0]


def read_lines(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return stream.readlines()


def write_lines(path, lines):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        stream.writelines(lines)


@dataclass(frozen=True)
class Element:
    n: int
    line: int                 # Zero-based input line; n is one-based.
    kind: str
    bpm: bool = False
    start: bool = False


@dataclass(frozen=True)
class Position:
    n: int
    kind: str
    z1: float
    z2: float


@dataclass(frozen=True)
class Pair:
    steerer: Element
    bpm: Element
    z: float


@dataclass(frozen=True)
class Monitor:
    element: Element
    z: float


@dataclass(frozen=True)
class Stage:
    label: int
    lines: tuple
    source: object           # None for START; SAVE label for RESTART.
    active: bool


def numeric_command(line):
    """Recognize commented numeric controls, excluding parameter legends."""
    text = line.lstrip().lstrip("!").lstrip()
    parts = code_part(text).split()
    required = {"SCHEFF": 7, "START": 5, "RESTART": 4, "SAVE": 0}
    if not parts or parts[0].upper() not in required:
        return None
    try:
        for value in parts[1:]:
            number(value)
    except ValueError:
        return None
    return text if len(parts) - 1 >= required[parts[0].upper()] else None


def stages(lines):
    """Catalog SCHEFF/tracking/SAVE triplets, including blocks after END.

    Uses the same catalog/activation approach as autophase2.stages(), retaining
    source line numbers so only genuine control commands are toggled.
    """
    result, pending, before_end = [], [], True
    for i, line in enumerate(lines):
        parts = code_part(line).split()
        if parts and parts[0].upper() in {"END", "END."}:
            before_end = False
        command = numeric_command(line)
        if command is None:
            if pending and parts:
                raise ValueError("Unexpected command inside tracking block at line %d" % (i + 1))
            continue
        words = code_part(command).split()
        key = words[0].upper()
        if key == "SCHEFF":
            if pending:
                raise ValueError("Incomplete tracking block before line %d" % (i + 1))
            pending = [i]
        elif key in {"START", "RESTART"}:
            if len(pending) != 1:
                raise ValueError("Each tracking stage requires SCHEFF before START/RESTART")
            pending.append(i)
        else:
            if len(pending) != 2:
                raise ValueError("Each SAVE requires one SCHEFF and one tracking command")
            label = number(words[1]) if len(words) > 1 else 0
            if label < 0 or label != int(label):
                raise ValueError("SAVE labels must be nonnegative integers")
            label = int(label)
            if any(s.label == label for s in result):
                raise ValueError("Duplicate SAVE label %d" % label)
            track = code_part(numeric_command(lines[pending[1]])).split()
            source = None
            if track[0].upper() == "RESTART":
                value = number(track[6]) if len(track) > 6 else 0
                if value < 0 or value != int(value):
                    raise ValueError("RESTART source must be a nonnegative SAVE label")
                source = int(value)
            active = before_end and bool(code_part(lines[pending[1]]).strip())
            result.append(Stage(label, tuple(pending + [i]), source, active))
            pending = []
    if pending:
        raise ValueError("Tracking block has no SAVE")
    if not result:
        raise ValueError("No SCHEFF / START-or-RESTART / SAVE stages found")
    return result


def select_stages(lines, catalog, first=None, last=None):
    """Activate an inclusive stage window and place END after its final SAVE.

    first=None selects geometry only: END precedes all tracking controls.
    Intermediate SAVE commands stay active because the following RESTART needs
    the newly tracked distribution, including the current trial correction.
    """
    if first is not None:
        if not 0 <= first <= last < len(catalog):
            raise ValueError("Invalid tracking stage window")
        writes = {s.label for s in catalog[first:last + 1]}
        if catalog[first].source in writes:
            raise ValueError("Stage window would overwrite its own upstream restart checkpoint")
        for index in range(first + 1, last + 1):
            if catalog[index].source != catalog[index - 1].label:
                raise ValueError("Stage SAVE %d must RESTART from preceding SAVE %d"
                                 % (catalog[index].label, catalog[index - 1].label))
    controls = {line_no: index for index, stage in enumerate(catalog) for line_no in stage.lines}
    insert_at = catalog[0].lines[0] if first is None else catalog[last].lines[-1] + 1
    result = []
    for i, line in enumerate(lines):
        if i == insert_at:
            result.append("END\n")
        if i in controls:
            command = numeric_command(line)
            if not command.endswith(("\n", "\r")):
                command += "\n"
            enabled = first is not None and first <= controls[i] <= last
            result.append(command if enabled else "!" + command)
        else:
            parts = code_part(line).split()
            result.append("!" + line if parts and parts[0].upper() in {"END", "END."} else line)
    if insert_at == len(lines):
        result.append("END\n")
    return result


@dataclass(frozen=True)
class Orbit:
    x: float
    y: float

    @property
    def score(self):
        return self.x * self.x + self.y * self.y

    @property
    def good(self):
        return abs(self.x) < BPM_TOLERANCE_MM and abs(self.y) < BPM_TOLERANCE_MM


def input_elements(lines, extra_names=()):
    """Count physical elements, including zero-length steerers and cathode."""
    names = ELEMENT_NAMES | set(extra_names)
    result = []
    title_payload = False
    for i, line in enumerate(lines):
        if title_payload:
            title_payload = False
            continue
        parts = code_part(line).split()
        if not parts:
            continue
        key = parts[0].lower()
        if key == "title":
            title_payload = True
            continue
        if key in {"end", "end."}:
            break
        if key not in names:
            continue
        if len(parts) < 2:
            raise ValueError("Incomplete %s element at input line %d" % (key, i + 1))
        number(parts[1])
        comment = line.partition(";")[2]
        start = bool(re.search(r"\bsteerer_start\b", comment, re.I))
        if start and key != "steerer":
            raise ValueError("steerer_start must mark a steerer (line %d)" % (i + 1))
        result.append(Element(len(result) + 1, i, key,
                              bool(re.fullmatch(r"\s*bpm\s*", comment, re.I)), start))
    return result


def select_pairs(elements):
    """Pair each selected steerer with the first subsequent BPM in the lattice."""
    starts = [e for e in elements if e.start]
    if len(starts) != 1:
        raise ValueError("Expected exactly one active 'steerer ... ; steerer_start' line; found %d"
                         % len(starts))
    bpms = [e for e in elements if e.bpm]
    used, pairs, skipped = set(), [], []
    for steerer in elements:
        if steerer.kind != "steerer" or steerer.n < starts[0].n:
            continue
        bpm = next((b for b in bpms if b.n > steerer.n), None)
        if bpm is None:
            skipped.append(steerer)
            continue
        if bpm.n in used:
            raise ValueError("Multiple selected steerers precede BPM at input line %d; "
                             "one-steerer/one-BPM pairing is ambiguous" % (bpm.line + 1))
        used.add(bpm.n)
        pairs.append((steerer, bpm))
    if not pairs:
        raise ValueError("No downstream '; bpm' element after steerer_start")
    return pairs, skipped


def select_svd_region(elements):
    """Independent corrector/BPM lists; SVD also supports rectangular systems."""
    starts = [e for e in elements if e.start]
    if len(starts) != 1:
        raise ValueError("Expected exactly one active 'steerer ... ; steerer_start' line; found %d"
                         % len(starts))
    monitors = [e for e in elements if e.bpm and e.n > starts[0].n]
    if not monitors:
        raise ValueError("No downstream '; bpm' element after steerer_start")
    candidates = [e for e in elements if e.kind == "steerer" and e.n >= starts[0].n]
    selected = [e for e in candidates if e.n < monitors[-1].n]
    skipped = [e for e in candidates if e.n >= monitors[-1].n]
    return selected, monitors, skipped


def read_positions(path):
    """Read the contiguous n/z1/element/z2 table, not later tracking tables.

    Steerer names can print as eight NUL bytes or blanks. Such rows still count;
    validation permits unnamed rows only for zero-length steerers.
    """
    rows, active = [], False
    with Path(path).open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if not active:
                if ELEMENT_HEADER.match(line):
                    active = True
                continue
            match = ELEMENT_ROW.match(line)
            if not match:
                if not rows and not line.strip():
                    continue
                break
            n, z1, name, z2, _ = match.groups()
            if int(n) != len(rows) + 1:
                raise ValueError("OUTPAR element indices are not consecutive at n=%s" % n)
            kind = name.replace("\x00", "").strip().lower()
            p = Position(int(n), kind, number(z1), number(z2))
            if p.z2 < p.z1:
                raise ValueError("OUTPAR element %s has z2 < z1" % n)
            rows.append(p)
    if not rows:
        raise ValueError("No n/z1/element/z2 table in %s" % path)
    return rows


def validate_elements(lines, positions):
    """Require exact element count and matching names in input and OUTPAR."""
    elements = input_elements(lines, {p.kind for p in positions if p.kind})
    if len(elements) != len(positions):
        raise ValueError("Element count mismatch: input=%d, OUTPAR=%d; "
                         "use OUTPAR from this exact input" % (len(elements), len(positions)))
    for e, p in zip(elements, positions):
        unnamed_steerer = not p.kind and e.kind == "steerer" and p.z1 == p.z2
        if e.kind != p.kind and not unnamed_steerer:
            raise ValueError("Element %d mismatch: input line %d is %s, OUTPAR is %s"
                             % (e.n, e.line + 1, e.kind, p.kind or "unnamed"))
    return elements


def map_elements(lines, positions):
    elements = validate_elements(lines, positions)
    raw_pairs, skipped = select_pairs(elements)
    pairs = []
    for steerer, bpm in raw_pairs:
        z = positions[bpm.n - 1].z2
        if z <= positions[steerer.n - 1].z2:
            raise ValueError("BPM n=%d is not longitudinally downstream of steerer n=%d"
                             % (bpm.n, steerer.n))
        pairs.append(Pair(steerer, bpm, z))
    return elements, pairs, skipped


def read_trajectory(path):
    """Read numeric samples after DATA; keep the last sample at duplicate Z."""
    columns, titles, active = [], False, False
    indices = None
    z_values, orbits = [], []
    with Path(path).open(encoding="utf-8", errors="replace") as stream:
        for line_no, line in enumerate(stream, 1):
            text = line.strip()
            if not active:
                if text.upper() == "TITLES":
                    titles, columns = True, []
                elif text.upper() == "ENDTITLES":
                    titles = False
                elif titles:
                    columns.append(re.sub(r"\s+", "", text).lower())
                elif text.upper() == "DATA":
                    active = True
                    try:
                        indices = tuple(columns.index(c) for c in
                                        ("z(cm)", "<x>(mm)", "<y>(mm)"))
                    except ValueError:
                        raise ValueError("Missing Z(cm), <X>(mm), or <Y>(mm) in table TITLES")
                continue
            if not text or text.startswith((";", "!", "#")):
                continue
            if text.upper() in {"END", "ENDDATA"}:
                break
            parts = text.split()
            try:
                z, x, y = (number(parts[i]) for i in indices)
            except (IndexError, ValueError) as exc:
                raise ValueError("Invalid trajectory data at line %d: %s" % (line_no, exc)) from exc
            if z_values and z < z_values[-1]:
                raise ValueError("Trajectory Z decreases at line %d; cannot mix overlapping "
                                 "tracking segments for interpolation" % line_no)
            if z_values and z == z_values[-1]:
                orbits[-1] = Orbit(x, y)
            else:
                z_values.append(z)
                orbits.append(Orbit(x, y))
    if not active or not z_values:
        raise ValueError("No numeric trajectory after DATA in %s" % path)
    return z_values, orbits


def orbit_at(trajectory, z):
    """Interpolate between the samples bracketing the BPM's z2."""
    zs, values = trajectory
    i = bisect_left(zs, z)
    if i < len(zs) and zs[i] == z:
        return values[i]
    if i == 0 or i == len(zs):
        raise ValueError("BPM Z=%.9g cm is outside table coverage [%.9g, %.9g] cm; "
                         "tracking must cover the BPM on both sides" % (z, zs[0], zs[-1]))
    weight = (z - zs[i - 1]) / (zs[i] - zs[i - 1])
    a, b = values[i - 1], values[i]
    return Orbit(a.x + weight * (b.x - a.x), a.y + weight * (b.y - a.y))


def steerer_values(line):
    parts = code_part(line).split()
    if len(parts) < 6 or parts[0].lower() != "steerer":
        raise ValueError("Invalid steerer line: %s" % line.rstrip())
    return number(parts[4]), number(parts[5])


def corrected_lines(lines, settings):
    """Change only two steerer values; preserve comments and all controls."""
    result = list(lines)
    for line_no, values in settings.items():
        line = result[line_no]
        spans = list(re.finditer(r"\S+", code_part(line)))
        steerer_values(line)
        for column, value in reversed(list(zip((4, 5), values))):
            if not math.isfinite(value):
                raise ValueError("Non-finite steerer amplitude")
            span = spans[column]
            line = line[:span.start()] + format(value, ".16g") + line[span.end():]
        result[line_no] = line
    return result


class Runner:
    def __init__(self, input_path, lines, executable="parmela", timeout=None, mode="1vs1"):
        self.input = Path(input_path).resolve()
        self.lines = lines
        self.executable, self.timeout = executable, timeout
        self.mode = mode
        method_suffix = "_svd" if mode == "svd" else "_1v1"
        self.trial = self.input.with_name(self.input.stem + "_temp.inp")
        self.output = self.input.with_name(self.input.stem + "_corrected" + method_suffix + ".inp")
        self.log_path = self.input.with_name(self.input.stem + "_autocorrection.log")
        self.summary = self.input.with_name(self.input.stem + "_correction_summary" + method_suffix + ".csv")
        self.outpar = self.input.parent / "OUTPAR.TXT"
        self.table = self.input.parent / "TIMESTEPEMITTANCE.TBL"
        self.cache, self.pairs, self.reference = {}, None, None
        self.count, self.checkpoint_written = 0, False
        self.catalog = stages(lines)
        self.first_stage = next((i for i, s in enumerate(self.catalog) if s.active), None)
        if self.first_stage is None:
            raise ValueError("No active tracking stage before END; enable the intended starting block")
        self.bounds = {}
        self.frontier = self.first_stage - 1
        self.window = None
        self.positions = None
        self.steerers, self.monitors, self.probe_stages = None, None, None
        self.last_full_svd = None
        self.log_path.write_text("PARMELA orbit correction (%s)\n" % mode, encoding="utf-8")

    def log(self, message):
        print(message, flush=True)
        with self.log_path.open("a", encoding="utf-8") as stream:
            stream.write(message + "\n")

    def checkpoint(self, settings):
        lines = corrected_lines(self.lines, settings)
        if self.window is not None:
            lines = select_stages(lines, self.catalog, *self.window)
        write_lines(self.output, lines)
        self.checkpoint_written = True

    def launch(self, lines, description, need_table=True):
        # Even a failed or shortened trial may overwrite SAVE distributions.
        # An orbit-cache hit alone can never establish physical SAVE validity.
        self.last_full_svd = None
        self.count += 1
        write_lines(self.trial, lines)
        # Never consume stale diagnostics if PARMELA fails to regenerate them.
        for path in (self.outpar, self.table):
            if path.exists():
                path.unlink()
        kwargs = {}
        if sys.platform == "win32":
            kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
        self.log("PARMELA run %d started at %s: %s" % (self.count, timestamp(), description))
        run_started = time.perf_counter()
        status = "failed"
        try:
            with self.log_path.open("a", encoding="utf-8") as stream:
                # A basename works for native Windows PARMELA from Cygwin Python too.
                subprocess.run([self.executable, self.trial.name], cwd=str(self.input.parent),
                               stdout=stream, stderr=subprocess.STDOUT, check=True,
                               timeout=self.timeout, **kwargs)
            status = "finished"
        except KeyboardInterrupt:
            status = "interrupted"
            raise
        finally:
            duration = time.perf_counter() - run_started
            self.log("PARMELA run %d %s at %s; elapsed %s (%.3f seconds)"
                     % (self.count, status, timestamp(), elapsed_text(duration), duration))
        if not self.outpar.is_file() or (need_table and not self.table.is_file()):
            raise ValueError("PARMELA did not produce fresh " +
                             ("OUTPAR.TXT and TIMESTEPEMITTANCE.TBL" if need_table else "OUTPAR.TXT"))
        positions = read_positions(self.outpar)
        if self.mode == "svd":
            elements = validate_elements(self.lines, positions)
            steerers, bpm_elements, skipped = select_svd_region(elements)
            monitors = [Monitor(e, positions[e.n - 1].z2) for e in bpm_elements]
            if monitors[0].z <= positions[steerers[0].n - 1].z2:
                raise ValueError("First SVD BPM must be longitudinally downstream of steerer_start")
            pairs = []
        else:
            elements, pairs, skipped = map_elements(self.lines, positions)
        geometry = [(p.kind, p.z1, p.z2) for p in positions]
        if self.reference is None:
            self.reference, self.pairs = geometry, pairs
            self.positions = positions
            if self.mode == "svd":
                self.steerers, self.monitors = steerers, monitors
                self.log("Validated %d input/OUTPAR elements; SVD uses %d steerers and %d BPMs."
                         % (len(elements), len(steerers), len(monitors)))
            else:
                self.log("Validated %d input/OUTPAR elements; %d steerer/BPM pairs."
                         % (len(elements), len(pairs)))
            for e in skipped:
                self.log("Skipping steerer n=%d, line %d: no downstream BPM."
                         % (e.n, e.line + 1))
        elif geometry != self.reference:
            raise ValueError("OUTPAR geometry changed between steerer trials")
        return read_trajectory(self.table) if need_table else None

    def initialize(self):
        """Obtain the full element table without tracking to the linac end."""
        self.launch(select_stages(self.lines, self.catalog),
                    "geometry only (END before all tracking blocks)", need_table=False)

    def execute_window(self, settings, first, last, description):
        lines = select_stages(corrected_lines(self.lines, settings), self.catalog, first, last)
        return self.launch(lines, description)

    def ensure_stage(self, index, settings):
        """Prepare upstream checkpoints once and learn actual Z coverage.

        Each preparation call tracks just one block. This avoids assuming SAVE
        labels are RF section numbers or estimating low-energy Z from phase.
        """
        while self.frontier < index:
            current = self.frontier + 1
            stage = self.catalog[current]
            if current > self.first_stage and stage.source != self.catalog[current - 1].label:
                raise ValueError("Stage SAVE %d must RESTART from preceding SAVE %d"
                                 % (stage.label, self.catalog[current - 1].label))
            trajectory = self.execute_window(settings, current, current,
                                             "prepare SAVE %d and measure stage coverage" % stage.label)
            self.bounds[current] = (trajectory[0][0], trajectory[0][-1])
            self.frontier = current
            self.log("  SAVE %d stage covers Z=%.9g to %.9g cm" %
                     (stage.label, *self.bounds[current]))

    def prepare_pair(self, pair, settings):
        """Find a restart before the steerer and include >=2 blocks when available."""
        steerer_z = self.positions[pair.steerer.n - 1].z2
        first = self.window[0] if self.window is not None else self.first_stage
        while True:
            self.ensure_stage(first, settings)
            if self.bounds[first][1] >= steerer_z:
                break
            first += 1
            if first == len(self.catalog):
                raise ValueError("Tracking stages never reach steerer Z=%.9g cm" % steerer_z)
        if self.bounds[first][0] > steerer_z:
            # A preceding block gives a safe upstream restart near a boundary.
            if first > self.first_stage:
                first -= 1
            else:
                raise ValueError("The first active stage starts downstream of steerer Z=%.9g cm; "
                                 "enable an earlier tracking block" % steerer_z)
        last = first
        while self.bounds[last][1] < pair.z or last == first:
            if last + 1 == len(self.catalog):
                if self.bounds[last][1] < pair.z:
                    raise ValueError("Tracking stages never reach BPM Z=%.9g cm" % pair.z)
                if first > self.first_stage:
                    first -= 1
                break
            last += 1
            self.ensure_stage(last, settings)
        self.window = (first, last)
        self.cache.clear()  # The restart checkpoint and target define a new objective.
        labels = [s.label for s in self.catalog[first:last + 1]]
        self.log("  Optimization window: SAVE blocks %s; END after SAVE %d."
                 % (", ".join(map(str, labels)), labels[-1]))
        if len(labels) == 1:
            self.log("  Only one tracking block is available; it covers this complete pair.")

    def run(self, settings, pair, force=False):
        key = tuple(sorted(settings.items()))
        if not force and key in self.cache:
            return self.cache[key]
        first, last = self.window
        trajectory = self.execute_window(settings, first, last,
                                         "%s SAVE %d through SAVE %d for BPM n=%d" %
                                         ("accepted settings / SAVE" if force else "trial",
                                          self.catalog[first].label, self.catalog[last].label,
                                          pair.bpm.n))
        result = orbit_at(trajectory, pair.z)
        self.cache[key] = result
        return result

    def prepare_svd(self, settings):
        # A common upstream anchor must precede EVERY variable corrector.
        self.prepare_pair(Pair(self.steerers[0], self.monitors[-1].element,
                               self.monitors[-1].z), settings)
        first, last = self.window
        self.probe_stages = []
        for steerer in self.steerers:
            z = self.positions[steerer.n - 1].z2
            index = first
            while index < last and self.bounds[index][1] < z:
                index += 1
            if self.bounds[index][0] > z and index > first:
                index -= 1
            self.probe_stages.append(index)
        if self.probe_stages != sorted(self.probe_stages):
            raise ValueError("SVD probe windows must follow longitudinal stage order")
        self.log("SVD global window: SAVE %d through SAVE %d. Response probes run "
                 "from downstream correctors toward upstream correctors, reusing baseline checkpoints."
                 % (self.catalog[first].label, self.catalog[last].label))

    def run_all(self, settings, force=False):
        key = tuple(sorted(settings.items()))
        if not force and key in self.cache:
            return self.cache[key]
        before = self.svd_output_versions()
        trajectory = self.execute_window(settings, *self.window,
                                         "SVD accepted settings / SAVE" if force else "SVD orbit measurement")
        result = {m.element.n: orbit_at(trajectory, m.z) for m in self.monitors}
        self.cache[key] = result
        manifest = self.svd_save_manifest()
        if manifest is not None:
            after = {name: version for name, version, _ in manifest}
            if before is None or any(after[name] == version for name, version in before.items()):
                manifest = None  # Every required output must have been regenerated.
        self.last_full_svd = (key, self.window, self.count, manifest, result)
        return result

    def svd_output_versions(self):
        if self.window is None:
            return None
        first, last = self.window
        versions = {}
        for stage in self.catalog[first:last+1]:
            if not 0 <= stage.label <= 26:
                return None
            name = "SAVECOR" + (chr(ord("a")+stage.label-1) if stage.label else "")
            try:
                stat = (self.input.parent / name).stat()
                versions[name] = (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
            except OSError:
                versions[name] = None
        return versions

    def svd_save_manifest(self):
        """Fingerprint the selected outputs AND their immutable input checkpoint.

        PARMELA's SAVE 0 is SAVECOR, SAVE 1 is SAVECORa, etc. Unknown naming or
        missing files disables the skip optimization instead of assuming safety.
        """
        if self.window is None:
            return None
        first, last = self.window
        labels = {s.label for s in self.catalog[first:last+1]}
        if self.catalog[first].source is not None:
            labels.add(self.catalog[first].source)
        manifest = []
        try:
            for label in sorted(labels):
                if not 0 <= label <= 26:
                    return None
                name = "SAVECOR" + (chr(ord("a")+label-1) if label else "")
                path = self.input.parent / name
                before = path.stat()
                digest = hashlib.sha256()
                with path.open("rb") as stream:
                    for block in iter(lambda: stream.read(1024*1024), b""):
                        digest.update(block)
                after = path.stat()
                version = lambda s: (s.st_size, s.st_mtime_ns, s.st_ctime_ns)
                if version(before) != version(after):
                    return None
                manifest.append((name, version(after), digest.hexdigest()))
        except OSError:
            return None
        return tuple(manifest)

    def ensure_svd_saved(self, settings, reason):
        """Skip a duplicate only if the LAST physical run saved this exact state."""
        saved = self.last_full_svd
        key = tuple(sorted(settings.items()))
        if saved is not None and saved[0] == key and saved[1] == self.window:
            if saved[3] is not None and saved[3] == self.svd_save_manifest():
                self.log("SVD %s: reusing verified full-window run %d; settings, window and "
                         "SAVE fingerprints match. Duplicate PARMELA run skipped."
                         % (reason, saved[2]))
                return saved[4]
        self.log("SVD %s: restoring accepted settings and SAVE files with a full-window run."
                 % reason)
        return self.run_all(settings, force=True)

    def probe_svd(self, settings, column, baseline):
        """Probe one amplitude; upstream BPMs retain their baseline orbit.

        Columns MUST be visited downstream-to-upstream. Each probe's source
        checkpoint then remains untouched by all preceding probe windows.
        """
        index, axis = divmod(column, 2)
        steerer = self.steerers[index]
        first = self.probe_stages[index]
        trajectory = self.execute_window(
            settings, first, self.window[1],
            "SVD response probe: steerer n=%d %s, SAVE %d through SAVE %d"
            % (steerer.n, "X" if axis == 0 else "Y", self.catalog[first].label,
               self.catalog[self.window[1]].label))
        result = list(baseline)
        for row, monitor in enumerate(self.monitors):
            if monitor.element.n > steerer.n:
                orbit = orbit_at(trajectory, monitor.z)
                result[2*row:2*row+2] = (orbit.x, orbit.y)
        return tuple(result)


def parabolic_axis(evaluate, current, axis, step):
    """Minimize x^2+y^2 with three samples and a safeguarded parabola vertex.

    evaluate((sx, sy)) returns this pair's interpolated BPM orbit.
    Only measured settings can be accepted, never an untested prediction.
    """
    center = current[axis]
    base = evaluate(current)
    if base.good:
        return current, base, step
    samples = [(current, base)]
    for direction in (-1, 1):
        trial = list(current)
        trial[axis] = center + direction * step
        if not math.isfinite(trial[axis]) or trial[axis] == center:
            raise ValueError("Scan step is unrepresentable at steerer amplitude %.9g" % center)
        trial = tuple(trial)
        orbit = evaluate(trial)
        samples.append((trial, orbit))
        if orbit.good:
            return trial, orbit, step
    f0, fm, fp = (samples[i][1].score for i in (0, 1, 2))
    if not all(math.isfinite(f) for f in (f0, fm, fp)):
        raise ValueError("Squared BPM orbit overflowed during parabolic fitting")
    # Fit in units of step, avoiding ill-conditioning for extreme amplitudes.
    a, b = (fm + fp) / 2 - f0, (fp - fm) / 2
    curvature_floor = 1e-12 * max(f0, fm, fp, 1e-300)
    if a > curvature_floor:
        offset = max(-MAX_STEP_MULTIPLIER,
                     min(MAX_STEP_MULTIPLIER, -b / (2 * a)))
        trial = list(current)
        trial[axis] = center + offset * step
        trial = tuple(trial)
        if math.isfinite(trial[axis]) and all(
                abs(trial[axis] - p[axis]) > 1e-12 * max(1, abs(center), step)
                for p, _ in samples):
            orbit = evaluate(trial)
            samples.append((trial, orbit))
            if orbit.good:
                return trial, orbit, step
    best, orbit = min(samples, key=lambda sample: sample[1].score)
    if orbit.score >= base.score:
        step *= 0.5
    elif abs(best[axis] - center) >= 0.95 * step:
        step *= 2
    return best, orbit, step


def optimize_pair(evaluate, initial, delta, iterations, log):
    """At most iterations X/Y sweeps, checking both planes after every trial."""
    current, steps = initial, [delta, delta]
    orbit = evaluate(current)
    if orbit.good:
        log("  Already within tolerance: X=%.6g, Y=%.6g mm; keeping amplitudes." % (orbit.x, orbit.y))
        return current, orbit, 0
    for iteration in range(1, iterations + 1):
        for axis in (0, 1):
            current, orbit, steps[axis] = parabolic_axis(evaluate, current, axis, steps[axis])
            log("  Sweep %d %s: steerer=(%.9g, %.9g), BPM=(%.6g, %.6g) mm"
                % (iteration, "X" if axis == 0 else "Y", *current, orbit.x, orbit.y))
            if orbit.good:
                return current, orbit, iteration
    log("  Iteration limit reached; retaining the best measured settings.")
    return current, orbit, iterations


def svd_correction(response, orbit, rcond=SVD_RCOND):
    """Return -R^+ orbit using a truncated one-sided Jacobi SVD.

    Orthogonalize R's columns while accumulating the right singular vectors V.
    The orthogonal columns are U*S. This avoids R.T*R (which squares the
    condition number) and needs no optional numerical package for small systems.
    """
    rows = len(response)
    columns = len(response[0]) if rows else 0
    if not rows or not columns or len(orbit) != rows:
        raise ValueError("SVD response/target dimensions are empty or inconsistent")
    if any(len(row) != columns for row in response):
        raise ValueError("SVD response matrix has inconsistent row lengths")
    if not 0 <= rcond < 1 or not math.isfinite(rcond):
        raise ValueError("SVD relative cutoff must be finite and in [0, 1)")
    if not all(math.isfinite(value) for row in response for value in row) or not all(
            math.isfinite(value) for value in orbit):
        raise ValueError("Non-finite SVD response or orbit")
    scale = max(abs(value) for row in response for value in row)
    if scale == 0:
        return (0.,)*columns, (0.,)*columns, 0
    # Store columns for the Jacobi rotations. Scale against overflow/underflow.
    b = [[response[i][j]/scale for i in range(rows)] for j in range(columns)]
    v = [[float(i == j) for i in range(columns)] for j in range(columns)]
    tiny = (sys.float_info.epsilon * max(rows, columns))**2
    for _ in range(100):
        rotations = 0
        for p in range(columns - 1):
            for q in range(p + 1, columns):
                alpha = math.fsum(x*x for x in b[p])
                beta = math.fsum(x*x for x in b[q])
                if min(alpha, beta) <= tiny:
                    continue
                gamma = math.fsum(x*y for x, y in zip(b[p], b[q]))
                if abs(gamma) <= 1e-12 * math.sqrt(alpha) * math.sqrt(beta):
                    continue
                zeta = (beta - alpha)/(2*gamma)
                t = math.copysign(1., zeta)/(abs(zeta) + math.hypot(1., zeta))
                c = 1/math.hypot(1., t)
                s = c*t
                b[p], b[q] = ([c*x-s*y for x, y in zip(b[p], b[q])],
                              [s*x+c*y for x, y in zip(b[p], b[q])])
                v[p], v[q] = ([c*x-s*y for x, y in zip(v[p], v[q])],
                              [s*x+c*y for x, y in zip(v[p], v[q])])
                rotations += 1
        if not rotations:
            break
    else:
        raise ValueError("Jacobi SVD did not converge within 100 sweeps")
    norms = [math.sqrt(math.fsum(x*x for x in column)) for column in b]
    cutoff = max(rcond, sys.float_info.epsilon * max(rows, columns)) * max(norms)
    retained = [j for j, norm in enumerate(norms) if norm > cutoff]
    factors = {j: -math.fsum(x*y for x, y in zip(b[j], orbit))/(norms[j]**2)/scale
               for j in retained}
    correction = tuple(math.fsum(v[j][i]*factors[j] for j in retained)
                       for i in range(columns))
    singular = tuple(sorted((norm*scale for norm in norms), reverse=True))
    if not all(math.isfinite(value) for value in correction):
        raise ValueError("SVD correction overflowed; increase the singular-value cutoff")
    return correction, singular, len(retained)


def optimize_svd(evaluate, initial, delta, iterations, log, rcond=SVD_RCOND,
                 probe=None, refresh=None, accept=None):
    """Reuse a measured response while checking every physical update.

    Poor predictions or failed progress on a reused matrix trigger one fresh
    response scan at the accepted baseline, within the same update iteration.
    Reverse probe order preserves the required upstream restart checkpoints.
    """
    current = tuple(initial)
    orbit = tuple(evaluate(current))
    if not orbit or not all(math.isfinite(x) for x in orbit):
        raise ValueError("Invalid SVD baseline orbit")
    rows = len(orbit)

    def good(values):
        return max(abs(x) for x in values) < BPM_TOLERANCE_MM

    def score(values):
        return math.fsum(x*x for x in values)

    def rms(values):
        return math.sqrt(score(values)/len(values))

    def checked(values):
        values = tuple(values)
        if len(values) != rows or not all(math.isfinite(x) for x in values):
            raise ValueError("Invalid orbit during SVD measurement")
        return values

    def result(values, measured, used):
        if accept is not None:
            accept(values)
        return values, measured, used

    log("SVD baseline: RMS=%.6g mm, max abs orbit=%.6g mm" %
        (rms(orbit), max(abs(x) for x in orbit)))
    if good(orbit):
        log("All selected BPMs already satisfy the tolerance; no response scan needed.")
        return result(current, orbit, 0)
    trust = MAX_STEP_MULTIPLIER*delta
    response, builds = None, 0
    for iteration in range(1, iterations + 1):
        best, best_orbit, best_score = current, orbit, score(orbit)
        best_reusable = False
        while True:  # At most a reused attempt followed by one fresh scan.
            reused = response is not None
            if not reused:
                if refresh is not None:
                    orbit = checked(refresh(current))
                    if good(orbit):
                        return result(current, orbit, iteration - 1)
                    if best == current or score(orbit) < best_score:
                        best, best_orbit, best_score = current, orbit, score(orbit)
                builds += 1
                log("SVD response matrix build %d for update iteration %d." % (builds, iteration))
                response = [[0.]*len(current) for _ in orbit]
                for column in reversed(range(len(current))):
                    trial = list(current)
                    trial[column] += delta
                    step = trial[column] - current[column]
                    if not math.isfinite(trial[column]) or step == 0:
                        raise ValueError("SVD scan step is unrepresentable at amplitude %.9g" % current[column])
                    trial = tuple(trial)
                    measured = checked(probe(trial, column, orbit) if probe is not None else evaluate(trial))
                    measured_score = score(measured)
                    if measured_score < best_score:
                        best, best_orbit, best_score = trial, measured, measured_score
                        best_reusable = False
                    if good(measured):
                        log("SVD response probe itself reached tolerance; retaining that measured setting.")
                        return result(trial, measured, iteration)
                    for row, (value, base) in enumerate(zip(measured, orbit)):
                        response[row][column] = (value - base)/step
            else:
                log("SVD iteration %d: reusing response matrix from build %d." % (iteration, builds))
            baseline_score = score(orbit)
            change, singular, rank = svd_correction(response, orbit, rcond)
            if not reused:
                log("SVD iteration %d: response %dx%d, retained rank %d, relative cutoff %.6g\n"
                    "  Singular values: %s" % (iteration, len(orbit), len(current), rank, rcond,
                                              ", ".join("%.6g" % s for s in singular)))
            if not rank:
                log("No usable response modes; retaining the best measured settings.")
                return result(best, best_orbit, iteration)
            largest = max(abs(x) for x in change)
            if largest > trust:
                change = tuple(x*(trust/largest) for x in change)
            rebuild = False
            for attempt in range(SVD_BACKTRACK_STEPS):
                gain = 0.5**attempt
                trial = tuple(value + gain*dx for value, dx in zip(current, change))
                if trial == current:
                    break
                measured = checked(evaluate(trial))
                measured_score = score(measured)
                movement = tuple(x-y for x, y in zip(trial, current))
                predicted_change = tuple(math.fsum(x*dx for x, dx in zip(row, movement))
                                         for row in response)
                prediction_error = math.sqrt(score(
                    tuple(value-base-dp for value, base, dp in zip(measured, orbit, predicted_change))))
                prediction_scale = max(math.sqrt(score(predicted_change)),
                                       0.1*BPM_TOLERANCE_MM*math.sqrt(rows), 1e-300)
                relative_error = prediction_error/prediction_scale
                accurate = relative_error <= SVD_PREDICTION_TOLERANCE
                log("  SVD gain %.6g: RMS=%.6g mm, max abs orbit=%.6g mm; "
                    "prediction error=%.3g%% (%s matrix)" %
                    (gain, rms(measured), max(abs(x) for x in measured), 100*relative_error,
                     "reused" if reused else "fresh"))
                if good(measured):
                    return result(trial, measured, iteration)
                if measured_score < best_score:
                    best, best_orbit, best_score = trial, measured, measured_score
                    best_reusable = accurate
                improved = measured_score < baseline_score*(1 - 1e-12)
                if reused and (not accurate or not improved):
                    log("SVD reused matrix failed prediction/progress check; rebuilding at "
                        "the accepted baseline before continuing this iteration.")
                    response, rebuild, best_reusable = None, True, False
                    break
                if improved:
                    break
            if reused and not rebuild and best_score >= baseline_score*(1 - 1e-12):
                log("SVD reused matrix produced no effective update; rebuilding once.")
                response, rebuild = None, True
            if rebuild:
                continue
            break
        if best_score >= baseline_score*(1 - 1e-12):
            log("SVD made no measured improvement; stopping at the best accepted orbit.")
            return result(current, orbit, iteration)
        movement = max(abs(x-y) for x, y in zip(best, current))
        if movement >= 0.95*trust:
            trust *= 2
        current, orbit = best, best_orbit
        if not best_reusable:
            response = None
            log("SVD response will be rebuilt: accepted setting lacks an accurate model prediction.")
        if accept is not None:
            accept(current)
        log("SVD accepted iteration %d: RMS=%.6g mm, max abs orbit=%.6g mm" %
            (iteration, rms(orbit), max(abs(x) for x in orbit)))
    log("SVD iteration limit reached; retaining the best measured settings.")
    return result(current, orbit, iterations)


def run_svd_correction(runner, args):
    controls = runner.steerers
    initial = tuple(value for e in controls for value in steerer_values(runner.lines[e.line]))

    def settings_for(values):
        return {e.line: tuple(values[2*i:2*i+2]) for i, e in enumerate(controls)}

    settings = settings_for(initial)
    runner.prepare_svd(settings)
    runner.checkpoint(settings)

    def vector(orbits):
        return tuple(value for m in runner.monitors
                     for value in (orbits[m.element.n].x, orbits[m.element.n].y))

    def evaluate(values):
        return vector(runner.run_all(settings_for(values)))

    def refresh(values):
        return vector(runner.ensure_svd_saved(settings_for(values), "response baseline"))

    def probe(values, column, baseline):
        return runner.probe_svd(settings_for(values), column, baseline)

    best, _, used = optimize_svd(evaluate, initial, args.delta_val, args.iterations,
                                 runner.log, args.svd_rcond, probe=probe, refresh=refresh,
                                 accept=lambda values: runner.checkpoint(settings_for(values)))
    settings = settings_for(best)
    # Restore after rejected/cached/shortened trials, but avoid repeating an
    # already verified full-window run with matching SAVE distributions.
    final = runner.ensure_svd_saved(settings, "final SAVE")
    runner.checkpoint(settings)
    first, last = (runner.catalog[index].label for index in runner.window)
    with runner.summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("bpm_n", "bpm_line", "bpm_z_cm", "orbit_x_mm", "orbit_y_mm",
                         "within_tolerance", "iterations", "first_save", "last_save", "mode",
                         "steerer_n", "steerer_line", "steerer_x", "steerer_y"))
        for m in runner.monitors:
            orbit = final[m.element.n]
            # A reference association for the CSV; SVD itself remains global.
            associated = next(e for e in reversed(controls) if e.n < m.element.n)
            writer.writerow((m.element.n, m.element.line + 1, m.z, orbit.x, orbit.y,
                             orbit.good, used, first, last, "svd", associated.n,
                             associated.line + 1, *settings[associated.line]))
            runner.log("Final SVD BPM n=%d: X=%.6g, Y=%.6g mm [%s]" %
                       (m.element.n, orbit.x, orbit.y, "OK" if orbit.good else "NOT CONVERGED"))
    for e in controls:
        runner.log("Final SVD steerer n=%d (line %d): X=%.9g, Y=%.9g" %
                   (e.n, e.line + 1, *settings[e.line]))
    success = all(orbit.good for orbit in final.values())
    runner.log("%s after %d PARMELA runs.\nCorrected input: %s\nSummary: %s" %
               ("All SVD BPMs within tolerance" if success else "SVD correction incomplete",
                runner.count, runner.output, runner.summary))
    return 0 if success else 2


def run_correction(args):
    correction_started = time.perf_counter()
    lines = read_lines(args.input)
    # Validate markers and pairing before starting an expensive simulation.
    selector = select_svd_region if args.mode == "svd" else select_pairs
    selector(input_elements(lines))
    runner = Runner(args.input, lines, args.parmela, args.timeout, mode=args.mode)
    settings, sweeps, final, windows = {}, {}, {}, {}
    try:
        runner.log("Correction started at %s" % timestamp())
        runner.log("Correction mode: %s" % args.mode)
        runner.log("Tolerance: abs(X) and abs(Y) < %.6g mm; absolute scan step: %.6g"
                   % (BPM_TOLERANCE_MM, args.delta_val))
        runner.log("Staged tracking: fixed upstream restart and accepted-settings SAVE runs.")
        runner.initialize()
        if args.mode == "svd":
            return run_svd_correction(runner, args)
        for i, pair in enumerate(runner.pairs, 1):
            line_no = pair.steerer.line
            initial = settings.get(line_no, steerer_values(lines[line_no]))
            runner.log("Pair %d/%d: steerer n=%d (line %d) -> BPM n=%d (line %d), Z=%.9g cm"
                       % (i, len(runner.pairs), pair.steerer.n, line_no + 1,
                          pair.bpm.n, pair.bpm.line + 1, pair.z))
            runner.prepare_pair(pair, settings)
            settings[line_no] = initial
            runner.checkpoint(settings)

            def evaluate(values):
                trial = dict(settings)
                trial[line_no] = values
                return runner.run(trial, pair)

            best, orbit, used = optimize_pair(evaluate, initial, args.delta_val,
                                               args.iterations, runner.log)
            settings[line_no], sweeps[pair.bpm.n] = best, used
            # Trials can leave SAVE files holding a rejected candidate. Always
            # regenerate the whole window at the accepted settings before moving.
            final[pair.bpm.n] = runner.run(settings, pair, force=True)
            windows[pair.bpm.n] = tuple(runner.catalog[index].label for index in runner.window)
            runner.checkpoint(settings)
        with runner.summary.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("steerer_n", "steerer_line", "bpm_n", "bpm_line", "bpm_z_cm",
                             "steerer_x", "steerer_y", "orbit_x_mm", "orbit_y_mm",
                             "within_tolerance", "sweeps", "first_save", "last_save"))
            for pair in runner.pairs:
                orbit = final[pair.bpm.n]
                writer.writerow((pair.steerer.n, pair.steerer.line + 1, pair.bpm.n,
                                 pair.bpm.line + 1, pair.z, *settings[pair.steerer.line],
                                 orbit.x, orbit.y, orbit.good, sweeps[pair.bpm.n],
                                 *windows[pair.bpm.n]))
                runner.log("Final BPM n=%d: X=%.6g, Y=%.6g mm [%s]"
                           % (pair.bpm.n, orbit.x, orbit.y, "OK" if orbit.good else "NOT CONVERGED"))
        success = all(orbit.good for orbit in final.values())
        runner.log("%s after %d PARMELA runs.\nCorrected input: %s\nSummary: %s"
                   % ("All BPMs within tolerance" if success else "Correction incomplete",
                      runner.count, runner.output, runner.summary))
        return 0 if success else 2
    except (OSError, ValueError, subprocess.SubprocessError, KeyboardInterrupt) as exc:
        runner.log("Stopped: %s\nDiagnostics: %s" % (exc or "interrupted", runner.log_path))
        if runner.checkpoint_written:
            runner.log("Last accepted correction checkpoint: %s" % runner.output)
        return 1
    finally:
        duration = time.perf_counter() - correction_started
        runner.log("Correction ended at %s" % timestamp())
        runner.log("Total correction time: %s (%.3f seconds)" % (elapsed_text(duration), duration))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="PARMELA .inp file")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("-1vs1", dest="mode", action="store_const", const="1vs1",
                       help="sequential one-corrector/next-BPM correction (default)")
    modes.add_argument("-svd", dest="mode", action="store_const", const="svd",
                       help="coupled SVD correction of all selected BPMs")
    parser.set_defaults(mode="1vs1")
    parser.add_argument("delta_val", type=float, help="positive absolute steerer scan step")
    parser.add_argument("iterations", type=int, help="maximum sweeps per pair (-1vs1), or global update steps (-svd)")
    parser.add_argument("--svd-rcond", type=float, default=SVD_RCOND,
                        help="discard singular modes below this fraction of the strongest (default: %(default)g)")
    parser.add_argument("--parmela", default="parmela", help="PARMELA executable (default: parmela)")
    parser.add_argument("--timeout", type=float, help="optional timeout in seconds for each PARMELA run")
    args = parser.parse_args(argv)
    if args.input.suffix.lower() != ".inp" or not args.input.is_file():
        parser.error("input must be an existing .inp file")
    if not math.isfinite(args.delta_val) or args.delta_val <= 0:
        parser.error("delta_val must be finite and positive")
    if args.iterations < 1:
        parser.error("iterations must be at least 1")
    if not math.isfinite(args.svd_rcond) or not 0 <= args.svd_rcond < 1:
        parser.error("svd-rcond must be finite and in [0, 1)")
    if args.timeout is not None and (not math.isfinite(args.timeout) or args.timeout <= 0):
        parser.error("timeout must be finite and positive")
    try:
        return run_correction(args)
    except (OSError, ValueError) as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
