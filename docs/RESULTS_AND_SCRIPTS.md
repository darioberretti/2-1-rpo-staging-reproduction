# Results and scripts

The six runners follow the scientific order of the paper and supporting study.
Run Setup before a scientific section. Figure numbers below are fixed within
each runner and remain unchanged when a section is run alone.

MAT results are stored under `results/{reference,generated}/<result folder>/`.
Figures and CSV tables are written under `outputs/{reference,generated}/`
in study folders numbered 01–06.

## [Run 01: reference RPO](../run_01_reference_rpo.m)

Result folder: `reference_rpo`. Figure folder: `01_reference_rpo`.

| Section / contribution | Supporting-study location | MAT result | Figures |
| --- | --- | --- | --- |
| A1 branches and reference RPO | Tables 3.2–3.3; Figures 3.1–3.2 | `a1_reference_rpo.mat` | 01–02 |
| Planar and spatial linear stability | Figures 3.5–3.6 | `linear_stability.mat` | 03–04 |
| Poincare section | Figure 3.7 | `poincare_section.mat` | 05 |
| Vertical perturbations and phase response | Figures 3.8–3.10 | `vertical_instability.mat` | 06–08 |

The supporting-study locations refer to the thesis. Its Chapter 3
Earth-access material is included in Run 03.

## [Run 02: RPO to deep space](../run_02_rpo_to_deep_space.m)

Result folder: `rpo_to_deep_space`. Figure folder: `02_rpo_to_deep_space`.

| Section / contribution | MAT result | Figures |
| --- | --- | --- |
| Encounter map, Earth-departure construction, Venus/Mars target arcs and altitude comparison | `method_exploration.mat` | 01–04 |
| Venus patched guess, shooting nodes and optimized transfer in both frames | `venus_validation.mat` | 05–08 |
| Mars patched guess, shooting nodes and optimized transfer in both frames | `mars_validation.mat` | 09–12 |
| Sun-phase flexibility: delta-v before/after retuning and flight time after retuning | `flexibility.mat` | 13–15 |
| Directional reachability and interpolated velocity cloud | `reachability.mat` | 16–17 |

The encounter map, target-arc panel, altitude comparison and before-retuning
plot are explanatory additions to the manuscript results.

## [Run 03: Earth to RPO](../run_03_earth_to_rpo.m)

Result folder: `earth_to_rpo`. Figure folder: `03_earth_to_rpo`.

| Section / contribution | Supporting-study location | MAT result | Figures |
| --- | --- | --- | --- |
| Earth-access maps and two connection geometries | Chapter 3 | `earth_access.mat` | 01–04 |
| Patched guess, shooting nodes and optimized low-Earth-orbit transfer in rotating/inertial views | Chapter 6 | `earth_to_rpo_validation.mat` | 05–15 |
| Departure-altitude continuation | Chapter 6 | `altitude_continuation.mat` | 16 |

## [Run 04: Sun–Earth L2 to RPO](../run_04_sel2_to_rpo.m)

Result folder: `sel2_to_rpo`. Figure folder: `04_sel2_to_rpo`.
Supporting study: Chapter 6.

| Section / contribution | MAT result | Figures |
| --- | --- | --- |
| Reference halo construction | `sel2_reference_halo.mat` | 01 |
| Both unstable-direction branches | `sel2_unstable_directions.mat` | 02 |
| Encounter-database proxy frontier | `database_proxy_pareto.mat` | 03 |
| Patched guess, shooting nodes and optimized transfer in both frames | `sel2_to_rpo_validation.mat` | 04–07 |
| Sun-phase continuation: cost, flight time and elapsed time | `sel2_to_rpo_continuation.mat` | 08–09 |

## [Run 05: stationkeeping](../run_05_stationkeeping.m)

Result folder: `stationkeeping`. Figure folder: `05_stationkeeping`.

| Section / contribution | MAT result | Figures / tables |
| --- | --- | --- |
| 1000 one-year trials for each periodic solution | `stationkeeping_comparison.mat` | Figure 01; setup and annual-cost CSV tables |

The CSV tables reproduce the paper's stationkeeping comparison.
The annual-cost distribution is a supporting figure.

## [Run 06: comet-like mission](../run_06_comet_case_study.m)

Result folder: `comet_case_study`. Figure folder: `06_comet_case_study`.
Supporting study: Chapter 7.

| Section / contribution | MAT result | Figures |
| --- | --- | --- |
| Warning-time law and probability tables | `warning_law.mat` | 01 |
| Combined L1/L2, L1, L2 and Earth accessibility: full six-year range and 400/800/1000-day views | `heliocentric_accessibility.mat` | 02–09 |
| RPO escape calibration and 37-direction boundary | `rpo_calibration.mat` | 10–11 |
| Shared-catalogue Monte Carlo architecture comparisons | `monte_carlo_analyses.mat` | 12–13 |

See [Reproducibility](REPRODUCIBILITY.md) for figure interpretation
and differences expected from fresh computation.
