# Project handoff

This mostly static repository reproduces the journal paper and supporting
study on a Moon-grazing Earth-Moon retrograde periodic orbit as a staging node.

Do not refactor scientific logic until its result provenance is recorded. MATLAB scripts must begin with `clc; close all; clear;` followed by `addpath(genpath('functions'));`; keep user settings at the top and code clearly commented.

All runners use ordered output folders under `outputs/reference/` and
`outputs/generated/`: `01_reference_rpo`, `02_rpo_to_deep_space`,
`03_earth_to_rpo`, `04_sel2_to_rpo`, `05_stationkeeping`, and
`06_comet_case_study`. Figure filenames begin with a fixed two-digit number
in scientific order, such as `01_a1_family_samples.png`. Define each new
runner's complete figure order in `rpo.figure_numbers` and use the shared
exporters; running a section independently must preserve its figure numbers.
Create runner figures with `rpo.create_figure`, supplying export dimensions
instead of setting figure `Position`. This keeps display tabs in MATLAB's
figure container and export sizing independent of the window size.
