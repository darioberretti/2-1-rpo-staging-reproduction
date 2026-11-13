# 2:1 RPO staging paper reproduction

MATLAB scripts accompanying Dario Berretti and Kenta Oshima’s paper,
*A Moon-grazing retrograde periodic orbit as a rapid-response deep-space
staging node: Assessment using an adaptable trajectory-design framework*.
The scripts recreate figures and tables from the paper and supporting study
and can rerun the numerical experiments for the 2:1 Earth–Moon retrograde
periodic orbit (RPO).

## Quick start

1. Download [`accepted_results.zip`](https://zenodo.org/records/22901811/files/accepted_results.zip?download=1)
   from [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.22901811).
2. Extract it into the repository root, where the `run_0*.m` files are.
   The ZIP already contains `results/reference/`; avoid an extra
   `accepted_results/` folder. If you already have a `results/` folder,
   you can extract elsewhere and copy only the extracted `reference/` folder
   into `results/`. Leave `results/generated/` untouched. If replacing an
   existing reference set, move the old `reference/` folder aside first.
3. Set MATLAB’s current folder to the repository root and run a script below.
   All six default to the Zenodo results
   (`workflow = "use_saved_results"`, `saved_results_source = "zenodo"`).
   To run an individual scientific section, run Setup first.

Replay needs only base MATLAB and the saved-results ZIP. It does not need
the database ZIPs or solvers. Figures and tables go to `outputs/reference/`.

| Runner | Study | Additional requirements for recomputation |
| --- | --- | --- |
| [Run 01](run_01_reference_rpo.m) | A1 family, reference RPO and stability | None |
| [Run 02](run_02_rpo_to_deep_space.m) | Deep-space transfers, flexibility and reachability | CasADi/IPOPT for validation and flexibility; Optimization Toolbox for reachability |
| [Run 03](run_03_earth_to_rpo.m) | Earth access and optimized low-Earth-orbit transfers | Optimization Toolbox for transfer solves |
| [Run 04](run_04_sel2_to_rpo.m) | Sun–Earth L2 halo access to the RPO | Optimization Toolbox for transfer solves |
| [Run 05](run_05_stationkeeping.m) | One-year stationkeeping Monte Carlo comparison | None |
| [Run 06](run_06_comet_case_study.m) | Comet-Interceptor-like mission study | Optimization Toolbox for RPO calibration |

See [Results and scripts](docs/RESULTS_AND_SCRIPTS.md) for the
section-to-contribution and figure mapping.

## Recompute or reuse results

Set `workflow = "run_experiment"` to recompute. Scientific settings are
frozen in [the configuration](config/reproduction_config.m).
Run 01 also needs reference results for its comparison stage.
For experiments that need a reference database, download the relevant ZIPs
listed in [Data](docs/DATA.md) and extract them into the repository root.
The Lambert solver is bundled. For Run 02 validation and flexibility,
install the CasADi MATLAB interface and add it to the MATLAB path.
Parallel Computing Toolbox is optional for supported sweeps and builders.

`save_computed_results = true` retains local MAT results.
To replay them, select `workflow = "use_saved_results"` and
`saved_results_source = "local_run"`.
Replacing existing generated results requires confirmation; batch execution
requires `overwrite_generated_results = true`. Run 05 always requires
that explicit setting and never prompts.

Figures are written to `outputs/{reference,generated}/` in numbered study
folders, such as `01_reference_rpo/`. Filenames have fixed two-digit
prefixes even when sections run independently. Choose PNG, PDF or both
with `figure_formats`. Retained MAT results go to `results/generated/`;
downloaded reference files are never overwritten.

Fresh studies can take minutes to hours. See
[Reproducibility](docs/REPRODUCIBILITY.md) for numerical differences
and interpretation limits.

## Rebuild databases

The [database builders](data_builders/) recreate the forward RPO–Moon,
Moon-to-RPO, and Sun–Earth L2-to-Moon encounter tables. They save to
`data/generated/` and can resume compatible checkpoints. Run 06 rebuilds
its heliocentric database in its accessibility section. See
[Data](docs/DATA.md) for source settings and builder instructions.
