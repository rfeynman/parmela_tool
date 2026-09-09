# autophase2.py

Use the original rr10.inp directly. Python 3.8+ and PARMELA are required.

Maximum energy (crest phase):

```sh
python autophase2.py rr10.inp --stage 4 --section 0 --mode energy
```

Minimum energy spread:

```sh
python autophase2.py rr10.inp --stage 4 --section 0 --mode spread
```

List available stage and RF section numbers:

```sh
python autophase2.py rr10.inp --list
```

Optimize all remaining linacs sequentially through RF section 13:

```sh
python autophase2.py rr10.inp --stage 4 --section 0 --mode spread --global
```

Use `--mode energy` for maximum-energy tuning instead. Both stage and section
increase by one: stage 4/section 0 through stage 17/section 13 in this example.
Starting at stage 5/section 1 processes sections 1 through 13. The script checks
the complete stage/checkpoint chain before launching any simulations.

Global mode automatically saves each section's optimized checkpoint before
starting the next section; `--save-checkpoint` is not needed. Trial phases still
suppress SAVE. The next stage uses the newly optimized incoming beam, and its
crest is measured again. All previously tuned phases remain in the lattice.

This is one sequential pass minimizing spread (or maximizing energy) at each
stage's endpoint, not a joint optimization against the final linac exit. A single
`--z` cannot describe multiple endpoints, so omit it for a multi-section global
run. The initial restart checkpoint must already exist and precede the first
selected RF block.

The same `rr10_optimized.inp` is updated after every successful section and
`rr10_summary.txt` records each section's crest, phase offset and measurements.
At completion, the optimized input contains every tuned phase and only the last
stage's tracking block is active. If a later stage fails, earlier completed
phases and checkpoints remain saved; the summary records the failure. There are
no per-section folders or extra result files. `rr10.inp` remains unchanged.

- `--stage 4` selects the block ending in SAVE 4, which restarts from SAVE 3.
- `--section 0` tunes the first CELL + TRWAVE RF block (zero-based, 0 through 13).
- These numbers are independent. Select an RF block that the restarted beam traverses.
- `--mode energy` finds the maximum output energy and stops at that crest.
- `--mode spread` minimizes spread first, then measures crest for the off-crest report and reruns the selected spread optimum.

## Automatic input selection and output files

The script reads all SCHEFF/START/RESTART/SAVE blocks, including commented blocks
and those after END. In the working copy it activates the selected stage,
comments the other tracking commands with `!`, comments old END commands, and
inserts END immediately after the selected SAVE. The full lattice and stage
catalog remain available. No manual comment editing or rr10v.inp is needed.

`rr10.inp` is never changed. Trials use `rr10_temp.inp`; after a successful final
run, this becomes `rr10_optimized.inp`. The only other optimizer output is
`rr10_summary.txt`, which records trial measurements, the crest phase, optimized
phase, signed off-crest angle, final energy/spread, and completion status.
In energy mode the off-crest angle is zero. In spread mode it is optimized phase
minus crest phase, wrapped to [-180,180) degrees, in the TRWAVE input convention.
Intentional phase differences within the RF block are preserved.

No timestamped folders, CSV, JSON, or separate crest/spread input files are created.
Repeated runs replace the previous optimized input and summary. The original
input remains the preserved source, so no redundant rr10_backup.inp is required.
If a run fails, the temporary input and failure summary remain for diagnosis;
an existing optimized input is replaced only after successful final evaluation.

All filenames derive from the input stem. Run from the simulation directory with
its field maps and restart checkpoints available. PARMELA regenerates its normal
output files in that directory; do not run concurrent scans there.

Trials suppress SAVE, preserving the restart checkpoints. The final evaluation
also suppresses SAVE unless `--save-checkpoint` is supplied. The resulting
optimized input always has the selected SAVE enabled: running it manually writes
that checkpoint. The original .inp file is still unchanged.

## Search controls

Both searches use local neighbors instead of a full-period scan. The default
initial phase step is 2 degrees: evaluate the starting phase and its +/-2-degree
neighbors, follow the improving direction, expand the step to bracket the nearby
optimum, then use bounded Brent refinement (parabolic interpolation with golden-section fallback). `--phase-step` changes the
initial step. The crest search has a maximum travel of +/-180 degrees around the
input phase; this is a bound, not a scan interval. Spread search defaults to the same one-period bounds around the original input phase and starts there when accelerating. Refinement tolerance is 0.05 degrees with at
most 40 iterations for walking and 40 for each refinement.

```sh
python autophase2.py rr10.inp --stage 4 --section 0 --mode spread --offset-range -30 30 --phase-tol 0.1
```

Spread candidates must always have output energy GREATER than incoming energy.
A candidate on the decelerating branch is rejected; the search refines back
toward feasible accelerating phases. The optional `--min-energy-fraction` now
defaults to 0 (no additional fraction-of-crest constraint). Set it to 0.9 if you
also want 90% of crest output energy. Even zero cannot disable the positive
energy-gain requirement. `--relative-spread` minimizes Del-kE/kE instead of the
default absolute Del-kE.

Before each section search, the script estimates incoming energy with a probe
from the same restart checkpoint: one integration step of 0.000001 RF degree,
output every step, SAVE disabled. This avoids assuming the ordinary stage table
contains entrance data. The estimate is just after the inlet, not an exact read
of the binary checkpoint. If PARMELA does not produce a table for this very short
probe, the script stops rather than guessing. For a single-section run you can
provide the known inlet energy instead, for example `--input-energy 61.3` (MeV).
Global mode measures each inlet independently. The summary records inlet energy
and final energy gain as well as crest and off-crest phase.

The parser reads kE(MeV) and Del-kE(MeV) by column name and uses the last valid
numeric table row. This is the end of recorded tracking, not an automatically
located physical linac exit. Use `--z 1572 --z-tolerance 0.5` to select the nearest
recorded row within 0.5 cm of Z=1572 cm. No interpolation is performed.

`--crest-span 60` restricts crest search to +/-60 degrees around the input phase.
If the initial phase is decelerating, the energy search walks toward increasing
energy before starting spread optimization. Invalid output or failed PARMELA
execution stops the search. A flat energy response also stops it:
check that the checkpoint precedes the RF block and tracking reaches its exit.
This deliberately finds a neighboring optimum, not a global optimum over all
RF phases. Inspect particle losses. Maximizing output energy corresponds to maximum
energy gain for a fixed incoming beam state.

See `python autophase2.py --help` for executable, timeout, and table options.

Under Cygwin the script runs `parmela rr10_temp.inp` from the input directory.
It passes the basename because native Windows PARMELA cannot read Cygwin
absolute paths such as `/home/.../rr10_temp.inp`. It waits for each run to finish;
the shell background operator `&` is not needed inside the optimizer.

## Manual reference

Parmela3.doc, chapter V section 38 (printed page 73), describes SAVE/RESTART as
the way to restart before a changed section without rerunning the whole beamline.
The manual does not document an IF/GOTO/skip-block command. Python handles the
comment selection; the resulting input uses ordinary PARMELA syntax.
Section 39's CONTINUE continues an existing state without SAVE/RESTART disk I/O;
it does not load a previous checkpoint or skip earlier tracking blocks.

Brent refinement reuses the measured bracketing points. Energy mode minimizes
negative output energy. Spread mode minimizes spread squared (or relative
spread squared with --relative-spread), preserving the same optimum and the
positive energy-gain constraint. Summary values remain spread, not its square.

In spread mode, the console labels MINIMUM ENERGY SPREAD, then the later
MAXIMUM ENERGY crest measurement used only for off-crest reporting. The final
run restores the selected spread phase. Explicit --offset-range bounds are
crest-relative and therefore require measuring crest first; an explicit positive
--min-energy-fraction or a decelerating starting phase also requires crest first.
Without these conditions, there is no preliminary energy optimization.
