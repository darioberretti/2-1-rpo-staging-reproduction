# Reproducibility

## Reference and generated results

Reference files are the accepted saved results. Replaying them recreates
figures without rerunning the experiments. Generated files are local
recomputations and never replace reference files.
The [README](../README.md) explains the workflow settings.

The runners use frozen settings in `config/reproduction_config.m`;
database builders keep their frozen recipe at the top of each builder.
Each study retains its own model normalization and units.
Fresh numerical solutions and runtimes can vary with floating-point
arithmetic and optimizer convergence.

## Interpretation and fresh-run differences

### Run 01: reference orbit

Stability figures show the largest nontrivial Floquet multiplier magnitude,
using the full complex multipliers. The Poincare section retains only
initial conditions with complete requested crossing chains.
The vertical phase figure reports position `z`.

### Run 02: deep-space transfers

Flexibility recomputation uses frozen optimization starts and a recorded
case-specific retuning sequence. This sequence repeats known solves;
it is not a general automatic tuner. Solver status is retained in the results.

Fresh directional reachability can converge to slightly different local optima.

The escape-direction comparison counts candidates passing a lunar-speed
range filter, not fully feasible flyby transfers.
Selecting a generated phase-aware database for analyses that used the
legacy dense table also changes their inputs; see [Data](DATA.md).

### Run 03: Earth access

The access maps are preliminary trajectory searches, not optimized
low-Earth-orbit transfers. Their C3 values are nondimensional.
The direct transfer solves at the configured departure altitude
(500 km by default); the separate sweep continues from 200 to 1000 km.

The saved representative transfer comes from the historical altitude
continuation. The fresh direct solve has a different initialization,
so its flight time, lunar minimum and shooting nodes can differ.
Reported lunar minima use adaptive trajectory samples, not continuous minimization.

### Run 04: Sun–Earth L2 access

The database frontier is a duration/excess-speed proxy, not an optimized
maneuver-cost frontier. Initial-guess plots propagate frozen unoptimized
shooting nodes and retain their discontinuities.

Fresh transfer solves can differ from saved reference solutions.
Exact agreement of the complete continuation with the current solver
has not been established. Waiting time is the halo departure-phase offset;
total elapsed time is waiting time plus flight time.

### Run 05: stationkeeping

Annual cost sums realized impulses applied before day 365; the final coast
can finish later. Completion means successful targeting and reaching
the horizon, not continuous confinement within 500 km.
Targeting residuals describe the estimated commanded trajectory before
execution error. Mean, P99 and the empirical distribution use finite
completed trials; failed trials remain in the saved result and completion count.

### Run 06: comet-like mission

Rebuilt Lambert tables can differ numerically from reference tables.
Fresh RPO calibration and directional-boundary solves are regenerated
solutions; bit-identical agreement with saved results is not promised.

Monte Carlo uses the published scalar RPO capability limit of 4.08 km/s,
independently of recalibration. Capability limits describe excess speed,
not additional maneuver cost or a guarantee at every departure epoch.
Failed scenarios remain in the probability denominator.
Unreached epoch quantiles are displayed at the six-year mission cap.
