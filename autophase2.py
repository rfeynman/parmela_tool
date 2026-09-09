#!/usr/bin/env python3
"""Find a TRWAVE section's energy crest and minimum-spread offset.



Standard library only. Run --help or see autophase2_README.md.
Stage selection understands the SCHEFF / START-or-RESTART / SAVE blocks
in rr10.inp, including commented numeric commands and blocks after END.
"""
import argparse
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def tokens(line):
    return re.split(r"[\s,=]+", re.split(r"[!;:]", line, 1)[0].strip())


def number(s):
    return float(s.replace("D", "E").replace("d", "e"))


def numeric_command(line):
    """Recognize actual tracking commands, not commented parameter legends."""
    line = line.lstrip().lstrip("!").lstrip()
    t = tokens(line)
    if t[0].upper() not in {"SCHEFF", "START", "RESTART", "SAVE"}:
        return None
    try:
        for v in t[1:]:
            number(v)
    except ValueError:
        return None
    required = {"SCHEFF": 7, "START": 5, "RESTART": 4, "SAVE": 0}
    return line if len(t) - 1 >= required[t[0].upper()] else None


def stages(lines):
    """Deliberately support the simple rr10 block layout; reject ambiguity."""
    result, pending, start = [], [], None
    for i, line in enumerate(lines):
        cmd = numeric_command(line)
        if cmd is None:
            if pending and tokens(line)[0]:
                raise ValueError("Unexpected command inside tracking block at line %d" % (i + 1))
            continue
        key = tokens(cmd)[0].upper()
        if key == "SCHEFF":
            if pending:
                raise ValueError("Incomplete tracking block before line %d" % (i + 1))
            pending, start = [cmd], i
        elif key in {"START", "RESTART"}:
            if len(pending) != 1:
                raise ValueError("Stage selector requires SCHEFF before tracking command")
            pending.append(cmd)
        elif key == "SAVE":
            if len(pending) != 2:
                raise ValueError("Stage selector requires one tracking command before SAVE")
            t = tokens(cmd)
            label = int(number(t[1])) if len(t) > 1 else 0
            if any(s[0] == label for s in result):
                raise ValueError("Duplicate SAVE number in stage catalog")
            result.append((label, start, pending + [cmd]))
            pending = []
    if pending:
        raise ValueError("Tracking block has no SAVE")
    return result


def select_stage(lines, label):
    catalog = stages(lines)
    selected = next((s for s in catalog if s[0] == label), None)
    if selected is None:
        raise ValueError("No stage ending in SAVE %s" % label)
    # Keep the full catalog, activate one block, and move END after its SAVE.
    index = catalog.index(selected)
    stop = catalog[index + 1][1] if index + 1 < len(catalog) else len(lines)
    result = []
    for i, line in enumerate(lines):
        cmd = numeric_command(line)
        if cmd is not None:
            active = selected[1] <= i < stop
            result.append(cmd if active else "!" + cmd)
            if active and tokens(cmd)[0].upper() == "SAVE":
                result.append("END")
        elif tokens(line)[0].upper() == "END":
            result.append("!" + line)
        else:
            result.append(line)
    return result


def active_input(lines):
    out = []
    for line in lines:
        out.append(line)
        if tokens(line)[0].upper() == "END":
            break
    if not out or tokens(out[-1])[0].upper() != "END":
        out.append("END")
    return out


def sections(lines):
    found = []
    for i, line in enumerate(lines[:-1]):
        if tokens(line)[0].upper() != "CELL":
            continue
        j, wave = i + 1, []
        while j < len(lines):
            key = tokens(lines[j])[0].upper()
            if not key:  # allow blank lines and comments inside a block
                j += 1
                continue
            if key != "TRWAVE":
                break
            wave.append(j)
            j += 1
        if wave:
            found.append((i, wave))
    return found


def replace_phase(line, value):
    # Replace only the phase token, retaining field coefficients and comments.
    spans = list(re.finditer(r"[^\s,=]+", re.split(r"[!;:]", line, 1)[0]))
    m = spans[4]
    return line[:m.start()] + format(value, ".12g") + line[m.end():]


def phased(lines, section, phase, save=False):
    result = list(lines)
    cell, waves = sections(lines)[section]
    origin = number(tokens(lines[cell])[4]) - 90
    delta = phase - origin
    # Preserve any intentional phase offsets between cells.
    for i in [cell] + waves:
        result[i] = replace_phase(lines[i], number(tokens(lines[i])[4]) + delta)
    if not save:
        result = ["!" + line if tokens(line)[0].upper() == "SAVE" else line for line in result]
    return result


def read_table(path, target_z=None, z_tolerance=0.5):
    headings, in_titles, rows = [], False, []
    for line in path.read_text(errors="replace").splitlines():
        text = line.strip()
        if text.upper() == "TITLES":
            in_titles = True
            continue
        if text.upper() == "ENDTITLES":
            in_titles = False
            continue
        if in_titles:
            headings.append(text)
            continue
        if not headings:
            continue
        try:
            row = [number(v) for v in text.split()]
        except ValueError:
            continue
        if len(row) == len(headings):
            rows.append(row)
    required = ["Z(cm)", "kE(MeV)", "Del-kE(MeV)"]
    if not rows or any(h not in headings for h in required):
        raise ValueError("No valid data or required TITLES in %s" % path)
    iz, ie, isp = [headings.index(h) for h in required]
    row = rows[-1] if target_z is None else min(rows, key=lambda r: abs(r[iz] - target_z))
    if not all(math.isfinite(v) for v in row):
        raise ValueError("Non-finite values in selected table row")
    if target_z is not None and abs(row[iz] - target_z) > z_tolerance:
        raise ValueError("Tracking did not reach requested Z within tolerance")
    if row[ie] <= 0 or row[isp] < 0:
        raise ValueError("Invalid energy/spread in table")
    return {"z_cm": row[iz], "energy_MeV": row[ie], "spread_MeV": row[isp],
            "relative_spread": row[isp] / row[ie]}


def bounded_brent(evaluate, lo, hi, samples, tolerance, max_iters):
    """Bounded Brent refinement: safeguarded parabolic steps, golden fallback.

    Reuse bracketing evaluations, retain boundary optima, and never fit through
    infinite penalties from the acceleration constraint. Tolerance is the final
    bracket width in degrees, not a change in energy or spread.
    """
    known = [(x, f) for x, f in samples.items() if lo <= x <= hi]
    if not known:
        x = (lo + hi) / 2
        known = [(x, evaluate(x))]
    known.sort(key=lambda pair: pair[1])
    x, fx = known[0]
    if not math.isfinite(fx):
        raise ValueError("Brent refinement requires a feasible starting phase")
    w, fw = known[1] if len(known) > 1 else (x, fx)
    v, fv = known[2] if len(known) > 2 else (w, fw)
    a, b = lo, hi
    d = e = 0.0
    golden = (3 - math.sqrt(5)) / 2
    for _ in range(max_iters):
        if b - a <= tolerance:
            return x, True
        midpoint = (a + b) / 2
        minimum_step = max(tolerance / 4, 1e-14 * max(1, abs(x)))
        parabolic = False
        if abs(e) > minimum_step and all(math.isfinite(f) for f in (fx, fw, fv)):
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2 * (q - r)
            if q > 0:
                p = -p
            q = abs(q)
            previous_e = e
            e = d
            if q > 0 and abs(p) < abs(q * previous_e / 2) and q*(a-x) < p < q*(b-x):
                d = p / q
                u = x + d
                if u-a < 2*minimum_step or b-u < 2*minimum_step:
                    d = math.copysign(minimum_step, midpoint-x)
                parabolic = True
        if not parabolic:
            e = b-x if x < midpoint else a-x
            d = golden * e
        u = x + (d if abs(d) >= minimum_step else math.copysign(minimum_step, d))
        u = min(b, max(a, u))
        if u == x:
            return x, b-a <= tolerance
        fu = evaluate(u)
        if math.isnan(fu):
            raise ValueError("Non-finite objective in Brent refinement")
        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, fv, w, fw, x, fx = w, fw, x, fx, u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv, w, fw = w, fw, u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu
    return x, b-a <= tolerance

def write_input(path, lines):
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def local_minimize(evaluate, start, lo, hi, step, tolerance, max_iters):
    """Bracket the neighboring minimum by walking downhill, then refine it."""
    start = min(hi, max(lo, start))
    points = sorted(set((max(lo, start-step), start, min(hi, start+step))))
    values = {x: evaluate(x) for x in points}
    best = min(points, key=values.get)
    if best == start:
        return bounded_brent(evaluate, points[0], points[-1], values, tolerance, max_iters)
    direction = 1 if best > start else -1
    previous, current = start, best
    for _ in range(max_iters):
        step *= 1.9
        nxt = min(hi, max(lo, current + direction*step))
        if nxt == current:
            return current, True
        value = evaluate(nxt)
        if value >= values[current]:
            values[nxt] = value
            return bounded_brent(evaluate, min(previous, nxt), max(previous, nxt),
                            values, tolerance, max_iters)
        values[nxt] = value
        previous, current = current, nxt
    return current, False


class Runner:
    def __init__(self, args, lines):
        self.args, self.lines = args, lines
        self.cache, self.count = {}, 0
        self.table = args.input.parent / args.table
        self.trial = args.input.with_name(args.input.stem + "_temp.inp")
        self.output = args.input.with_name(args.input.stem + "_optimized.inp")
        self.summary = args.input.with_name(args.input.stem + "_summary.txt")
        if self.table.resolve() in {args.input, self.trial, self.output, self.summary}:
            raise ValueError("Table path conflicts with an input or summary file")
        self.summary.write_text(
            "PARMELA phase optimization\nAlgorithm: local bracketing + bounded Brent refinement\nStatus: RUNNING\n"
            "Input: %s\nMode: %s\nStage (destination SAVE): %d\nRF section: %d\n"
            "Crest search: input phase +/- %g deg\nSpread bounds: %g to %g deg (input-relative by default; crest-relative with --offset-range)\n"
            "Minimum fraction of crest output energy: %g\n"
            "Measurement: %s\nOriginal input is never modified.\n\n" %
            (args.input.name, args.mode, args.stage, args.section, args.crest_span,
             *(args.offset_range or (-args.crest_span, args.crest_span)), args.min_energy_fraction,
             "last valid table row" if args.z is None else "nearest Z=%g cm" % args.z))

    def log(self, text):
        with self.summary.open("a") as f:
            f.write(text + "\n")
        print(text, flush=True)

    def run(self, phase, force=False, save=False):
        key = round(phase, 10)
        if key in self.cache and not force:
            return self.cache[key]
        self.count += 1
        write_input(self.trial, phased(self.lines, self.args.section, phase, save))
        # PARMELA regenerates this output. Never accept a previous run's table.
        if self.table.exists():
            self.table.unlink()
        kwargs = {}
        if sys.platform == "win32":
            kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
        with tempfile.TemporaryFile(mode="w+") as stream:
            try:
                # PARMELA is a native Windows executable. Cygwin absolute paths
                # (/home/...) are not valid Windows filenames. Since cwd is the
                # input directory, a basename works in both Cygwin and Windows.
                # Wait for completion: shell '&' would race the table reader.
                subprocess.run([self.args.parmela, self.trial.name], cwd=str(self.args.input.parent),
                               stdout=stream, stderr=subprocess.STDOUT, check=True,
                               timeout=self.args.timeout, **kwargs)
            except (OSError, subprocess.SubprocessError):
                stream.seek(0)
                # Only retain diagnostics on failure, in the same summary file.
                self.log("PARMELA diagnostics:\n" + stream.read()[-8000:])
                raise
            if not self.table.exists():
                stream.seek(0)
                self.log("PARMELA diagnostics:\n" + stream.read()[-8000:])
                raise ValueError("PARMELA did not produce a new table")
        result = read_table(self.table, self.args.z, self.args.z_tolerance)
        self.cache[key] = result
        self.log("Run %d: phase=%.6f deg, E=%.9g MeV, spread=%.9g MeV, Z=%.6f cm" %
                 (self.count, phase, result["energy_MeV"], result["spread_MeV"], result["z_cm"]))
        return result


def optimize(args, lines, runner):
    origin = number(tokens(lines[sections(lines)[args.section][0]])[4])-90
    input_energy = args.input_energy
    if input_energy is None:
        # Probe just after the restart to measure the incoming beam using the
        # same table diagnostic. Ordinary stage tables may omit early rows.
        probe = list(lines)
        for i, line in enumerate(probe):
            t = tokens(line)
            if t[0].upper() not in {"START", "RESTART"}:
                continue
            offset = 1 if t[0].upper() == "RESTART" else 2
            t[offset], t[offset+1], t[offset+3] = "0.000001", "1", "1"
            probe[i] = " ".join(t)
        old_z = args.z
        try:
            runner.lines, args.z = probe, None
            runner.log("Measuring incoming energy with a one-step 0.000001-degree probe (SAVE disabled).")
            input_energy = runner.run(origin, force=True)["energy_MeV"]
        finally:
            runner.lines, args.z = lines, old_z
            runner.cache.clear()
    runner.log("Incoming energy: %.12g MeV" % input_energy)
    warnings = []
    crest = crest_data = None
    crest_ok = True
    def find_crest():
        runner.log("Search objective: MAXIMUM ENERGY (crest reference).")
        phase, ok = local_minimize(lambda x: -runner.run(x)["energy_MeV"], origin,
            origin-args.crest_span, origin+args.crest_span, args.phase_step, args.phase_tol, args.max_iters)
        data = runner.run(phase)
        if data["energy_MeV"] <= input_energy:
            raise ValueError("No accelerating crest found within the local search bounds")
        energies = [r["energy_MeV"] for r in runner.cache.values()]
        if max(energies)-min(energies) <= max(1e-10, abs(data["energy_MeV"])*1e-8):
            raise ValueError("Energy is flat versus phase. Check RF section, restart checkpoint and tracking endpoint.")
        return phase, data, ok

    if args.mode == "energy":
        crest, crest_data, crest_ok = find_crest()
        best, refined = crest, crest_ok
    else:
        # Only explicit crest-dependent constraints or a decelerating initial
        # state require crest-first work. Normally optimize spread immediately.
        initial = runner.run(origin)
        if args.min_energy_fraction > 0 or args.offset_range is not None or initial["energy_MeV"] <= input_energy:
            runner.log("Crest required first: explicit crest constraint or nonaccelerating starting phase.")
            crest, crest_data, crest_ok = find_crest()
        threshold = input_energy if crest_data is None else max(input_energy, crest_data["energy_MeV"] * args.min_energy_fraction)
        metric = "relative_spread" if args.relative_spread else "spread_MeV"
        def objective(x):
            r = runner.run(x)
            if r["energy_MeV"] <= input_energy or r["energy_MeV"] < threshold:
                runner.log("Reject phase %.6f: insufficient accelerating energy." % x)
                return math.inf
            return r[metric] ** 2
        lo, hi = (origin-args.crest_span, origin+args.crest_span) if args.offset_range is None else [crest+v for v in args.offset_range]
        start = min(hi, max(lo, origin))
        if not math.isfinite(objective(start)) and crest is not None:
            start = min(hi, max(lo, crest))
        if not math.isfinite(objective(start)):
            raise ValueError("No feasible accelerating starting phase in the requested bounds")
        runner.log("Search objective: MINIMUM ENERGY SPREAD (accelerating branch).")
        best, spread_ok = local_minimize(objective, start, lo, hi, args.phase_step, args.phase_tol, args.max_iters)
        if not math.isfinite(objective(best)):
            raise ValueError("No accelerating spread minimum found")
        runner.log("Minimum-spread phase selected: %.9f deg; spread=%.12g MeV." % (best, runner.run(best)["spread_MeV"]))
        if min(abs(best-lo), abs(best-hi)) <= args.phase_tol:
            warnings.append("Spread optimum is at a search boundary; consider widening the bounds.")
        if crest is None:
            runner.log("Measuring crest ONLY for off-crest reporting; minimum-spread phase remains selected.")
            crest, crest_data, crest_ok = find_crest()
        refined = spread_ok and crest_ok
    if not refined:
        warnings.append("Refinement iteration limit reached; increase --max-iters.")
    if args.crest_span < 180 and min(abs(crest-origin+args.crest_span), abs(crest-origin-args.crest_span)) <= args.phase_tol:
        warnings.append("Crest is at search boundary; increase --crest-span.")
    # A final run puts the requested optimum in PARMELA's normal output files.
    runner.log("Final run: restore selected %s optimum at %.9f deg." % (args.mode, best))
    final = runner.run(best, force=True, save=args.save_checkpoint)
    if final["energy_MeV"] <= input_energy:
        raise ValueError("Final verification is not accelerating; optimized input not published")
    # Publish only after successful final evaluation. Keep the original input intact.
    write_input(runner.trial, phased(lines, args.section, best, save=True))
    runner.trial.replace(runner.output)
    offset = (best-crest+180)%360-180
    runner.log(("\nStatus: SECTION COMPLETE\n" if args.global_run else "\nStatus: COMPLETE\n") +
        "Optimized input: %s\nCrest phase: %.9f deg\n"
        "Crest output energy: %.12g MeV\nOptimized phase: %.9f deg\n"
        "Off-crest phase (optimized minus crest): %+.9f deg\n"
        "Final output energy: %.12g MeV\nFinal energy spread: %.12g MeV\n"
        "Final relative energy spread: %.12g\nEvaluation Z: %.9f cm\n"
        "Refinement tolerance reached: %s\nCheckpoint saved during final run: %s\n"
        "The optimized input includes SAVE; running it manually updates that checkpoint." %
        (runner.output.name, crest, crest_data["energy_MeV"], best, offset,
         final["energy_MeV"], final["spread_MeV"], final["relative_spread"], final["z_cm"],
         refined, args.save_checkpoint))
    for warning in warnings:
        runner.log("Note: " + warning)
    runner.log("Energy gain: %.12g MeV" % (final["energy_MeV"] - input_energy))
    return phased(lines, args.section, best, save=True)


def global_schedule(lines, first_stage, first_section):
    """One stage per RF section, with a verified consecutive checkpoint chain."""
    if not 0 <= first_section <= 13 or len(sections(lines)) < 14:
        raise ValueError("--global requires 14 RF sections and a starting section from 0 to 13")
    catalog = {label: cmds for label, _, cmds in stages(lines)}
    schedule = []
    for section in range(first_section, 14):
        stage = first_stage + section - first_section
        if stage not in catalog:
            raise ValueError("Missing stage %d required for RF section %d" % (stage, section))
        track = tokens(catalog[stage][1])
        if track[0].upper() != "RESTART":
            raise ValueError("Global linac optimization requires a RESTART for stage %d" % stage)
        source = int(number(track[6])) if len(track) > 6 else 0
        if source != stage - 1:
            raise ValueError("Stage %d must RESTART from SAVE %d; found SAVE %d" % (stage, stage-1, source))
        schedule.append((stage, section))
    return schedule


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path)
    p.add_argument("--section", type=int, help="zero-based CELL + TRWAVE block; required for optimization")
    p.add_argument("--stage", type=int, help="destination SAVE number; automatically activates this tracking stage")
    p.add_argument("--mode", choices=("energy", "spread"), help="maximize energy or minimize spread; required for optimization")
    p.add_argument("--global", dest="global_run", action="store_true",
                   help="optimize successive stage/section pairs through section 13; saves each optimized checkpoint")
    p.add_argument("--list", action="store_true", help="list RF sections and tracking stages without simulation")
    p.add_argument("--crest-span", type=float, default=180, help="crest search +/- degrees around input phase (default 180)")
    p.add_argument("--offset-range", type=float, nargs=2, default=None, metavar=("LOW", "HIGH"), help="optional spread bounds relative to crest; supplying this requires crest-first search")
    p.add_argument("--phase-step", type=float, default=10, help="initial local neighbor step in degrees (default 2)")
    p.add_argument("--input-energy", type=float, help="known stage input energy in MeV; otherwise measured with a tiny tracking probe")
    p.add_argument("--phase-tol", type=float, default=0.05, help="refinement interval tolerance in degrees")
    p.add_argument("--max-iters", type=int, default=40, help="maximum refinement iterations per sampled basin")
    p.add_argument("--relative-spread", action="store_true", help="minimize Del-kE/kE instead of Del-kE")
    p.add_argument("--min-energy-fraction", type=float, default=0, help="optional fraction of crest energy; net acceleration is always required (default 0)")
    p.add_argument("--z", type=float, help="evaluate table row nearest this Z in cm; default last valid row")
    p.add_argument("--z-tolerance", type=float, default=0.5)
    p.add_argument("--parmela", default="parmela", help="PARMELA executable path")
    p.add_argument("--table", default="TIMESTEPEMITTANCE.TBL", help="table filename relative to simulation directory")
    p.add_argument("--timeout", type=float, default=1800, help="seconds allowed per simulation")
    p.add_argument("--save-checkpoint", action="store_true", help="execute SAVE on the final optimized run only")
    args = p.parse_args(argv)
    args.input = args.input.resolve()
    if args.input.suffix.lower() != ".inp":
        p.error("input must have an .inp extension")
    lines = args.input.read_text().splitlines()
    if args.list:
        for n, (i, wave) in enumerate(sections(lines)):
            print("RF section %d: line %d, %d TRWAVE lines, phase %.8f deg" %
                  (n, i+1, len(wave), number(tokens(lines[i])[4])-90))
        for label, _, cmds in stages(lines):
            print("Stage %d: %s -> %s" % (label, cmds[1], cmds[2]))
        return
    if args.stage is None or args.section is None or args.mode is None:
        p.error("optimization requires --stage, --section and --mode energy|spread")
    lines = select_stage(lines, args.stage)
    if not 0 <= args.section < len(sections(lines)):
        p.error("--section is outside the detected RF section range")
    values = [args.crest_span, args.phase_step, *(args.offset_range or ()), args.phase_tol, args.z_tolerance,
              args.timeout, args.min_energy_fraction]
    if not all(math.isfinite(x) for x in values) or (args.z is not None and not math.isfinite(args.z)):
        p.error("numeric settings must be finite")
    if not 0 < args.crest_span <= 180 or args.phase_step <= 0 or args.phase_tol <= 0 or args.max_iters < 1:
        p.error("require 0 < crest-span <= 180, phase-step > 0, phase-tol > 0, max-iters >= 1")
    if args.input_energy is not None and (not math.isfinite(args.input_energy) or args.input_energy <= 0):
        p.error("--input-energy must be positive and finite")
    if (args.offset_range is not None and args.offset_range[0] >= args.offset_range[1]) or not 0 <= args.min_energy_fraction <= 1 or args.z_tolerance < 0 or args.timeout <= 0:
        p.error("invalid offset bounds, minimum energy fraction, Z tolerance or timeout")
    reads = {int(number(tokens(l)[6])) if len(tokens(l)) > 6 else 0 for l in lines if tokens(l)[0].upper() == "RESTART"}
    writes = {int(number(tokens(l)[1])) if len(tokens(l)) > 1 else 0 for l in lines if tokens(l)[0].upper() == "SAVE"}
    if reads & writes:
        p.error("SAVE would overwrite the checkpoint used by RESTART")
    schedule = [(args.stage, args.section)]
    if args.global_run:
        schedule = global_schedule(lines, args.stage, args.section)
        if args.z is not None and len(schedule) > 1:
            p.error("--z is a single endpoint; omit it for --global so each stage uses its own final table row")
        if args.input_energy is not None and len(schedule) > 1:
            p.error("omit --input-energy in multi-section --global; each stage measures its own incoming energy")
        args.save_checkpoint = True
    runner = Runner(args, lines)
    try:
        if args.global_run:
            runner.log("Global sequence: " + ", ".join("stage %d / section %d" % pair for pair in schedule))
            runner.log("Each optimized stage is saved before the next stage starts.")
        for stage, section in schedule:
            args.stage, args.section = stage, section
            lines = select_stage(lines, stage)
            runner.lines = lines
            runner.cache.clear()  # A new checkpoint and RF block define a new objective.
            runner.log("\nOptimizing stage %d / RF section %d" % (stage, section))
            lines = optimize(args, lines, runner)
        if args.global_run:
            runner.log("\nStatus: COMPLETE\nGlobal sequence finished through RF section 13.\n"
                       "All optimized phases are retained in %s; its active tracking block is stage %d." %
                       (runner.output.name, args.stage))
    except (ValueError, OSError, subprocess.SubprocessError, KeyboardInterrupt) as exc:
        runner.log("\nStatus: FAILED\n%s\nTemporary input retained: %s" % (exc, runner.trial))
        raise


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print("autophase2: %s" % exc, file=sys.stderr)
        sys.exit(1)
