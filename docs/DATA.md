# Data

Download files from [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.22901811).
Extract each ZIP into the repository root, preserving the folder paths
inside it. The saved-results ZIP installs `results/reference/`; the
database ZIPs install `data/reference/`. Leave local computations in
`results/generated/` and `data/generated/` untouched.

| ZIP | Needed for |
| --- | --- |
| `accepted_results.zip` | Saved-result replay for all six runners; reference comparisons in Run 01 and reference Earth-access inputs in Run 03 |
| `rpo_to_moon_phase_aware.zip` | Run 02 diagnostics, validation and flexibility; Run 06 RPO calibration |
| `rpo_to_moon_legacy.zip` | Run 02 escape-direction comparison and reachability with the reference dense table |
| `sel2_to_moon.zip` and `moon_to_rpo.zip` | Run 04 database proxy |
| `comet_lambert_databases.zip` | Run 06 accessibility and Monte Carlo |

Only `accepted_results.zip` is needed for replay. Zenodo also provides
`manifest.json` and `SHA256SUMS` for verifying downloads.

## Saved results and small inputs

The [script map](RESULTS_AND_SCRIPTS.md) lists each section's MAT result.
Its path is `results/{reference,generated}/<study>/<filename>`, where
`<study>` is the result folder in that map.
Replay needs the selected saved result and the inputs already included
in the repository. No database download is needed for replay.

Small inputs are included in `data/inputs/`:

- `reference_rpo_seed.mat`: planar RPO state;
- `rpo_moon_database_seed.mat`: frozen RPO state shared by both encounter-table builders;
- `rpo_to_deep_space_targets.mat`: Venus/Mars target vectors;
- `run_02_flexibility_initial_guesses.mat` and
  `run_02_reachability_initial_guesses.mat`: frozen optimization starts;
- `run_03_reference_inputs.mat`, `run_04_reference_inputs.mat` and
  `run_05_reference_inputs.mat`: study-specific initial conditions and setups.

Frozen optimization starts contain unoptimized problem definitions.
The optimizations are still executed during recomputation.

## Forward RPO–Moon databases: Runs 02 and 06

Reference files under `data/reference/rpo_to_moon/`:

| File | Role |
| --- | --- |
| `RPOMoonDatabase10x3_50phases.mat` | Phase-aware table for Run 02 diagnostics, validation and flexibility, and Run 06 calibration |
| `RPOMoonDatabase50x3.mat` | Legacy dense table for Run 02 escape-direction comparison and reachability |

Run [the forward builder](../data_builders/build_rpo_to_moon_database.m)
to create the phase-aware table from the tracked seed. It needs no downloaded
database or earlier runner. Its frozen grid uses 50,000 departures,
0–50 m/s impulses and an 85-day horizon.
The table and metadata are saved under `data/generated/rpo_to_moon/`.

Select `experiment_database_source = "generated"` to use that table.
Run 02 then uses it for every section, including the two analyses that used
the legacy dense table in the reference study; differences are expected.
The builder does not regenerate the legacy dense table.

Compatible progress in `cache/rpo_to_moon_checkpoint.mat` resumes
automatically. If the builder reports incompatible progress, explicitly
remove that named checkpoint to restart. Existing generated data requires
one replacement confirmation, or `overwrite_generated_database = true`
in batch execution. The completed file is replaced only after a successful build.

The table uses 384405 km per distance unit and 4.34811305 days per time unit.
Lunar entry is at Moon radius plus Laplace sphere-of-influence radius,
about 67928 km from the Moon. `vInfMoon_RTH` is the local relative
velocity at that boundary, not an asymptotic hyperbolic velocity.

## Earth-access tables: Run 03

Fresh transfer and altitude-continuation sections load `earth_access.mat`
from `results/{reference,generated}/earth_to_rpo/`.
Choose the branch with `experiment_database_source`.
For generated inputs, first compute and save the Earth-access section.
No separate database builder is needed.

## Encounter tables: Run 04

Only the fresh database-proxy section needs these files:

- `data/reference/sel2_to_moon/SEL2MoonDatabase_BCR4BP_50000_1e-6_thetaSun0.mat`;
- `data/reference/moon_to_rpo/RPOMoonDatabase10x3_50phases_backward_5dv.mat`.

`experiment_database_source = "generated"` selects the same filenames
under `data/generated/`; both files must be present there.
Run [the Moon-to-RPO builder](../data_builders/build_moon_to_rpo_database.m)
to create the backward table from the tracked RPO seed: 50 phases, 0–5 m/s
impulses and 200 days of backward propagation. Event times remain negative;
the proxy uses their magnitudes. Its checkpoint is
`cache/moon_to_rpo_checkpoint.mat`, with the same resume and replacement
controls as the forward builder.

Run [the SEL2–Moon builder](../data_builders/build_sel2_to_moon_database.m)
to create the other table. It reconstructs the halo internally and uses the
tracked Run 04 direction to fix the unstable branch's sign; no earlier run
or downloaded database is needed. It propagates 50,000 plus-branch seeds
for up to 120 days in Sun–Earth CR3BP and 500 additional days in Earth–Moon
BCR4BP. The output contains all lunar entries and compact diagnostics,
without dense trajectory histories or plots. Compatible progress resumes
from `cache/sel2_to_moon_checkpoint.mat`, with the same replacement controls.

The SEL2 table uses 384400 km per Earth–Moon distance unit and
`(1/2.665e-6)/86400` days per Earth–Moon time unit. Seed Sun phase advances
with elapsed halo time. `t_enc_rel_days_em` includes both propagation stages
but excludes the initial halo-phase wait; absolute encounter times include
that wait. `vInfMoon_RTH_km_s` is the local relative velocity at the lunar
Laplace sphere-of-influence boundary. Run 04 applies its separate 500-day
cap to the sum of SEL2 and backward-RPO encounter durations.

Fresh transfer solves use the tracked
unoptimized shooting-node problems independently of this database choice.
Run 05 needs no database.

## Heliocentric databases: Run 06

Fresh accessibility and Monte Carlo sections need
`data/reference/comet_lambert/comet_lambert_database.mat`.
It contains the shared radius, angle and flight-time axes and cumulative
departure excess-speed grids for Earth, L1 and L2. The analyses derive the
combined L1/L2 envelope from those grids.
The database is unnecessary for saved-result replay.

Set `heliocentric_database_mode` to `"zenodo"` for reference tables,
`"generated"` for local tables, or `"rebuild"` to create them in the
accessibility section. Run that section before Monte Carlo when rebuilding.
The generated destination is `data/generated/comet_lambert/`;
replacement requires confirmation or `overwrite_generated_results = true`.

RPO calibration independently uses `experiment_database_source`.
Both controls default to `"generated"` for recomputation. To use the
downloaded reference inputs, set each relevant control to `"zenodo"`.
