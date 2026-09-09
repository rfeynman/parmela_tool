#!/usr/bin/env python3
"""Apply quadrupoles.csv column J to the 13 inter-structure PARMELA quads.

Usage: python update_quadrupoles.py INPUT.inp
Updates INPUT.inp in place, saving the original as INPUT.inp.bak (or .bak.N).
Searches beside the input, then beside this script, for quadrupoles.csv or
Sband14_output/quadrupoles.csv. Only Python's standard library is required.
"""

import argparse
import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path
import re


def find_csv(inp):
    for directory in (inp.parent, Path(__file__).resolve().parent):
        for relative in ("quadrupoles.csv", "Sband14_output/quadrupoles.csv"):
            candidate = directory / relative
            if candidate.is_file():
                return candidate
    raise ValueError("Cannot find quadrupoles.csv or Sband14_output/quadrupoles.csv "
                     "beside the input file or script.")


def read_settings(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = csv.reader(stream)
        header = next(rows, [])
        if len(header) < 10 or header[9].strip() != "Gradient_Gauss_per_cm":
            raise ValueError("CSV column J must be headed Gradient_Gauss_per_cm.")
        settings = []
        for line_number, row in enumerate(rows, 2):
            if not any(field.strip() for field in row):
                continue
            if len(row) < 10:
                raise ValueError(f"CSV row {line_number} has no column J.")
            value = row[9].strip()
            try:
                number = Decimal(value)
            except InvalidOperation:
                raise ValueError(f"Invalid gradient at CSV row {line_number}: {value!r}")
            if not number.is_finite():
                raise ValueError(f"Non-finite gradient at CSV row {line_number}.")
            if row[0].strip() != f"Q{len(settings) + 1}":
                raise ValueError("CSV rows must be ordered Q1 through Q13.")
            settings.append(str(number))
    if len(settings) != 13:
        raise ValueError(f"Expected 13 CSV settings; found {len(settings)}.")
    return settings


def replace_settings(data, settings):
    # Latin-1 provides a reversible byte mapping, preserving existing encoding.
    lines = data.decode("latin-1").splitlines(keepends=True)
    records = []
    for index, line in enumerate(lines):
        code = re.split(r"[!;]", line, maxsplit=1)[0]
        tokens = list(re.finditer(r"\S+", code))
        if tokens:
            records.append((index, tokens[0].group().lower(), tokens))

    # A structure starts with CELL followed by TRWAVE, ignoring blank/comments.
    starts = [i for i in range(len(records) - 1)
              if records[i][1] == "cell" and records[i + 1][1] == "trwave"]
    if len(starts) != 14:
        raise ValueError(f"Expected 14 CELL/TRWAVE structures; found {len(starts)}.")
    if len(settings) != 13:
        raise ValueError("Expected exactly 13 gradient settings.")

    for gap, (start, stop) in enumerate(zip(starts, starts[1:])):
        end = start + 1
        while end < stop and records[end][1] == "trwave":
            end += 1
        quads = [record for record in records[end:stop] if record[1] == "quad"]
        if len(quads) != 1:
            raise ValueError(f"Expected one quad between structures {gap + 1} and "
                             f"{gap + 2}; found {len(quads)}.")
        index, _, tokens = quads[0]
        if len(tokens) < 5:
            raise ValueError(f"Quad at input line {index + 1} has no amplitude.")
        amplitude = tokens[4]
        lines[index] = (lines[index][:amplitude.start()] + settings[gap]
                        + lines[index][amplitude.end():])
    return "".join(lines).encode("latin-1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inp", type=Path, help="PARMELA input file to update in place")
    args = parser.parse_args()
    try:
        inp = args.inp.resolve()
        csv_path = find_csv(inp)
        original = inp.read_bytes()
        updated = replace_settings(original, read_settings(csv_path))
        if updated == original:
            print(f"All 13 gradients already match {csv_path}; no changes needed.")
            return
        backup = inp.with_name(inp.name + ".bak")
        suffix = 0
        while True:
            try:
                with backup.open("xb") as stream:
                    stream.write(original)
                break
            except FileExistsError:
                suffix += 1
                backup = inp.with_name(inp.name + f".bak.{suffix}")
        inp.write_bytes(updated)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Error: {exc}\n")
    print(f"Updated 13 quad amplitudes (Gauss/cm) in {inp}\n"
          f"CSV: {csv_path}\nBackup: {backup}")


if __name__ == "__main__":
    main()
