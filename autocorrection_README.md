# PARMELA orbit correction: one-to-one and SVD

Run from Cygwin or another Python 3.7+ environment (standard library only):

```sh
python autocorrection.py rr11.inp -1vs1 0.1 10
python autocorrection.py rr11.inp -svd 0.1 10
```

Use `-1vs1` for the established sequential correction or `-svd` for coupled
correction of the selected BPMs together. Omitting the mode retains the previous
one-to-one behavior, so `python autocorrection.py rr11.inp 0.1 10` still works.

`delta_val` is a **positive absolute steerer amplitude scan step**. `iterations`
limits X/Y fitting sweeps **per pair** in `-1vs1` mode, or global correction
update iterations in `-svd` mode. A matrix rebuild, when needed, happens within
the same update iteration. Unlike
the old script, a nonzero initial amplitude does not make the step relative.
There is no `start_sect` argument or `!cor` requirement.

Optional arguments:

```sh
python autocorrection.py rr11.inp -1vs1 0.1 10 --parmela /path/to/parmela --timeout 3600
```

## Selecting the correction region and BPMs

Put exactly one `steerer_start` tag after the semicolon on an active steerer:

```text
steerer 0 4.5 1 0 0 ; steerer_start
```

That steerer and all later active steerers are candidates. Earlier steerers
remain unchanged. In `-1vs1`, each candidate controls the **first subsequent** element
whose semicolon comment is exactly `bpm` (case insensitive):

```text
drift 0.1 1.5 1 ; bpm
```

The measurement location is the **end** of that element, not its midpoint.
Comments such as `; bpmcube` do not create a monitor. A steerer without a
downstream BPM is logged and skipped. Multiple selected steerers sharing the
same next BPM are rejected because the requested one-to-one pairing is ambiguous.
In `-svd`, correctors and monitors are independent lists: every tagged BPM after
`steerer_start` is included, and all candidate correctors before the last BPM
are variables. Different corrector/BPM counts are supported. Correctors after
the last BPM are skipped in both modes.
Commented lattice elements and lattice elements after the first active `END`
are ignored. The tracking-block catalog also reads commented numeric controls
and later blocks after `END`, so the optimization window can advance.

## Element mapping and orbit measurement

An initial **geometry-only** PARMELA run puts `END` before all tracking blocks
and generates fresh `OUTPAR.TXT` without tracking to the linac end.
The script counts physical input elements in file order,
including zero-length `cathode` and `steerer` elements. Field setup commands
such as `poisson`, `cfield`, and field filenames do not occupy these indices.
Additional named element types are discovered from OUTPAR's element table.
The input's complete element count and sequence must match the contiguous
`n / z1 / element / z2` table. Blank/NUL element names are accepted only for
zero-length steerers, as printed by the PARMELA build used here.

Each BPM's `z2` supplies its longitudinal location in cm. The script reads
numeric timestep rows **after `DATA`**, identifies `Z(cm)`, `<X>(mm)`, and
`<Y>(mm)` through `TITLES`, and linearly interpolates between the two samples
bracketing the BPM. Exact sample positions use that sample directly. The last
data row is included. Duplicate Z values retain the last sample; decreasing Z,
malformed data, and BPMs outside the sampled range cause an explicit error.
The code does not extrapolate or substitute the nearest endpoint.

## One-to-one fitting and stopping

`BPM_TOLERANCE_MM = 0.03` near the top of `autocorrection.py` sets acceptance
for both modes (the configured value is preserved):

```text
abs(X at BPM) < 0.03 mm AND abs(Y at BPM) < 0.03 mm
```

This is a per-plane limit, not a radial-distance limit. A pair already within
the limit needs no scan. It is not necessary to reach exact zero.

For the current pair only, the objective is `F = X**2 + Y**2`. Each sweep:

1. Holds the Y amplitude fixed and samples F at X amplitude `u-h`, `u`, `u+h`.
2. Fits `F(u+h*t) = A*t**2 + B*t + C`; a positive-curvature fit predicts
   `t_min = -B/(2*A)`.
3. Caps the prediction at four current scan steps and simulates it. Only an
   actually measured improvement is retained. Flat/downward fits use the best
   sampled setting instead.
4. Repeats for the Y amplitude, then repeats X/Y sweeps if necessary. Both orbit
   components enter the objective so the method can handle moderate X/Y coupling.

The scan step grows when the best setting moves to or beyond a scan edge and
shrinks when there is no improvement. Both orbit limits are checked after each
measurement, so a satisfactory setting stops the pair immediately. Earlier
accepted steerers remain applied as correction moves downstream. The search is
local; it does not guarantee convergence for strong coupling, nonlinear response,
inadequate steerer response, or unattainable orbit targets. There are no physical
amplitude bounds in this input format; the four-step cap limits each prediction,
not the total accumulated amplitude.

## SVD correction

SVD uses a joint orbit vector `[X_BPM1, Y_BPM1, X_BPM2, Y_BPM2, ...]` and
amplitude vector `[X_steerer1, Y_steerer1, X_steerer2, Y_steerer2, ...]`. A
forward probe of `delta_val` measures each response-matrix column:

```text
R[:, j] = (orbit(amplitudes + delta_val * e_j) - baseline_orbit) / delta_val
R = U S V^T
amplitude_change = -V S_truncated_inverse U^T baseline_orbit
```

This includes X/Y cross-plane coupling and each corrector's effect on every
downstream BPM. Upstream BPM responses to downstream correctors are zero.
The truncated pseudoinverse minimizes the sum of squared BPM X/Y offsets;
for an underdetermined system it gives the minimum-norm retained-mode update.
Reaching the configured tolerance at every BPM is still required for success.

The production solver uses a standard-library one-sided Jacobi SVD, without
forming `R^T R` and without requiring NumPy or SciPy. `SVD_RCOND = 1e-3` discards
modes whose singular value is at or below 0.001 times the largest singular
value. A machine-precision floor also applies. To change the relative cutoff:

```sh
python autocorrection.py rr11.inp -svd 0.1 10 --svd-rcond 0.001
```

The response matrix, retained rank, singular values, gains, and measured orbit
statistics are described in the log (matrix dimensions are logged, not all
matrix entries). The initial maximum amplitude change is four scan steps. This
limit can grow after successful large updates. A predicted update is simulated;
if necessary up to six gains are tested, from 1 down to 1/32. Only a measured improvement
or a setting satisfying every BPM's tolerance is accepted. The measured matrix
is reused for later updates while its predictions remain accurate; incrementing
the iteration number does not by itself trigger another response scan.

For each simulated update the code compares the measured orbit change with
`R * applied_amplitude_change`. A relative prediction error above
`SVD_PREDICTION_TOLERANCE = 0.25` (25%), or lack of measured improvement with a
reused matrix, triggers a fresh scan at the accepted baseline during the same
iteration. The error denominator is the predicted change's Euclidean norm,
with a floor of `0.1 * BPM_TOLERANCE_MM * sqrt(number_of_orbit_components)` to
avoid treating negligible changes as large relative errors. A fresh matrix's
poor prediction also marks it for replacement after an improving step. A step
that actually reaches every BPM's tolerance succeeds regardless of prediction
error. Flat/ineffective responses can stop early; the iteration limit is an
upper bound. Matrix builds, reuse, prediction errors and rebuild reasons are
recorded in the log.

SVD tracking covers the first selected corrector through the last selected BPM,
starting from a common fixed upstream checkpoint. This is a larger region than
the one-to-one window. To reduce response-scan time, correctors are probed from
**downstream to upstream**. Each probe starts at the tracking block upstream
of its own corrector and tracks through the last selected BPM, while upstream
BPMs keep their baseline readings. This order prevents a probe from overwriting
the upstream baseline checkpoint needed by a subsequent probe. Before a new
response scan, the accepted baseline's SAVE state is checked and restored only
if the most recent physical run did not already produce that state.

With M correctors there are up to 2M response probes per matrix build, plus
baseline and damped-update measurements. An update using an accurate existing
matrix requires no response probes. SVD is not guaranteed to be faster than one-to-one
correction. There is no automatic full-cathode-to-linac-end scan: upstream
preparation and future block truncation follow the same stage catalog as before.

Accepted SVD amplitudes are checkpointed during optimization. Before baseline
restoration or final SAVE, the code checks whether the **last successful physical
full-window run** already used those exact settings and window, measured all
selected BPMs, and regenerated every required SAVE file. It also compares SHA-256
content hashes and file metadata for the output SAVE files and upstream source
checkpoint. If everything matches, it reuses that run and logs the skipped
duplicate. A cached orbit by itself is insufficient. Other launches, including
shortened probes and failed trials, invalidate this proof. Changed, missing,
unrefreshed or unrecognized SAVE files require a full-window restoration run.
Supported SAVE names are PARMELA's `SAVECOR`, `SAVECORa` through `SAVECORz`;
other naming disables skipping conservatively.

The original input is unchanged. The summary CSV has one row per BPM, with
final X/Y orbit, convergence, global iteration count, window labels, and
`mode=svd`. It also records `steerer_x` and `steerer_y`, plus `steerer_n` and
`steerer_line` to identify those settings. The association is the nearest
selected upstream corrector for reporting; SVD still corrects all BPMs globally.
One-to-one summaries already include these corrector settings and retain their
existing columns. All corrector amplitudes are also recorded in the log and
the corrected input. The output filenames,
per-run timestamps, total timing, and success/failure exit statuses are shared
with one-to-one mode.

## One-to-one tracking and common outputs

The stage catalog uses the same approach as `autophase2.py`: recognize actual
numeric `SCHEFF / START-or-RESTART / SAVE` triplets, including commented blocks
and blocks after `END`. Parameter legends and ordinary comments stay comments.
The first originally active tracking block is the starting point. If it is a
RESTART, its input saved distribution must already exist. Upstream blocks that
are still needed are run once to prepare the initial correction region; there
is no full-linac baseline run.

Preparation tracks one block at a time and records its actual timestep Z range.
This identifies which block contains the steerer and which reaches the BPM,
without assuming SAVE labels correspond to RF section numbers or converting
phase to distance at an assumed beam energy. These ranges are learned only as
far downstream as currently needed.

For each pair, the script selects a window whose restart is upstream of the
steerer and whose end reaches the BPM. It keeps **at least two consecutive
blocks when available**, extending the window further if needed to reach the
BPM. A file with only one usable block can use that block if it covers the
whole pair. Completed upstream and future `SCHEFF`, tracking, and `SAVE`
commands are commented, and `END` is inserted immediately after the last
selected `SAVE`. Lattice definitions are retained.

For example, a two-block window ending in SAVE 4 and SAVE 5 activates only:

```text
SCHEFF ...
RESTART ... 3
SAVE 4
SCHEFF ...
RESTART ... 4
SAVE 5
END
```

All trials restart from the same upstream SAVE 3 distribution. Intermediate
SAVE commands remain active: the second block must read the first block's
current trial distribution. The window is rejected if it would overwrite its
own upstream restart checkpoint or break the consecutive SAVE/RESTART chain.

After finding acceptable settings (or reaching the iteration limit), the script
**reruns that complete window at the accepted amplitudes**, verifies the BPM,
and writes its SAVE distributions. Only then does it advance. The next window
may overlap the last block because the next steerer can lie inside that block.
It reuses the accepted upstream checkpoint, rather than the last rejected trial.
The objective cache is cleared when the stage window/checkpoint changes.

There is **no full-linac final verification run**. The summary records each
pair's accepted-settings SAVE run. Subsequent downstream steerer changes do not
alter already corrected upstream BPMs in this forward-tracking setup.

`iterations` limits X/Y sweeps, not subprocess count. Each axis uses up to three
new window runs per sweep; cached repeats are reused. There are also one-time
stage preparation runs and one accepted-settings SAVE run per pair. Existing
generated `OUTPAR.TXT` and `TIMESTEPEMITTANCE.TBL` are removed before each launch
so failures cannot silently reuse stale diagnostics.

For `rr11.inp`, the outputs are:

- `rr11_corrected_svd.inp` or `rr11_corrected_1v1.inp`: accepted settings and the current stage window, with
  completed/future blocks commented and END after that window's final SAVE.
- `rr11_temp.inp`: the most recent PARMELA trial; after normal completion it
  matches the accepted final settings.
- `rr11_correction_summary_svd.csv` or `rr11_correction_summary_1v1.csv`: accepted SAVE-run verification of every paired
  BPM, element indices, source line numbers, z2, amplitudes, orbit, convergence,
  and the first/last SAVE labels in its optimization window.
- `rr11_autocorrection.log`: progress and PARMELA stdout/stderr diagnostics.

The corrected input and summary use `_svd` for `-svd` and `_1v1` for `-1vs1`
(also for the default mode when the flag is omitted), immediately before the
extension. This keeps results from the two methods separate. Temporary input
and log filenames retain their existing names.

The log records each PARMELA run's start/end timestamp (local time with UTC
offset) and elapsed duration. Its final line reports the total correction time
in `HH:MM:SS.mmm` and seconds, including geometry, preparation, optimization,
and accepted-settings SAVE runs. Timings are also recorded when a run fails or
is interrupted. Durations use a monotonic clock, so clock adjustments do not
affect elapsed-time measurements.

The original input is never edited. If a pair hits the iteration limit, the
best measured setting is retained, its window is rerun to SAVE that state, and
processing continues downstream. Exit status is 0
when every paired BPM passes, 2 when correction finishes with remaining orbit
errors, and 1 on simulation/parsing failure or interruption. Argument errors
also use argparse's exit status 2.

## Verification

```sh
python -m unittest -v test_autocorrection
```

The tests cover indexing, blank/NUL steerer names, strict pairing, interpolation,
reordered table columns, Fortran exponents, tolerance, coupled planes, iteration
limits, stage selection and END placement, stale-output rejection, and the
complete sequential workflow with a simulated PARMELA process and persistent
SAVE distributions. SVD tests cover coupled and rectangular systems, weak-mode
truncation, nonlinear damping, probe ordering, immutable restart checkpoints,
conditional matrix reuse and fallback, iteration limits, and final accepted
SAVE states. Tests also cover missing/tampered/unrefreshed SAVE files, cached
orbits after rejected trials, changed windows and failed launches. When NumPy is available the
solver is also compared against its pseudoinverse; this reference check is
optional and NumPy is not needed to run either correction mode. The tests do
not launch the real simulation.
