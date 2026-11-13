clc; close all; clear;
addpath(genpath('functions'));

%% Scientific purpose and instructions
% Reproduce Earth access and optimized low-Earth-orbit-to-RPO transfers,
% producing 16 figures.
%
% Run from the repository root. Run Setup before individual scientific sections.
% Saved-result replay uses three installed MAT results and base MATLAB.
% See docs/DATA.md for input locations.
%
% workflow = "run_experiment" recomputes with config/reproduction_config.m.
% The million-point access search supports optional parallel workers.
% Transfer solves require Optimization Toolbox and saved Earth-access tables;
% compute and save the access section first when selecting generated inputs.
% The direct transfer altitude defaults to 500 km; the separate continuation
% uses its own altitude grid.
%
% Figures: outputs/{reference,generated}/03_earth_to_rpo/.
% Retained MAT results: results/generated/earth_to_rpo/.
% See README.md for shared controls and docs/REPRODUCIBILITY.md for
% fresh-run differences.

%% User settings
workflow = "use_saved_results";        % "use_saved_results" or "run_experiment"

% Choose the saved result source for replay.
saved_results_source = "zenodo";       % "zenodo" or "local_run"

% Choose Earth-access inputs and result retention for fresh computation.
experiment_database_source = "zenodo"; % "zenodo" or "generated" Earth-access tables
save_computed_results = true;
overwrite_generated_results = false;

figure_formats = ["png"];              % "png", "pdf", or both

%% Setup and configuration

% Locate the repository and load the frozen scientific settings.
repository_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repository_root,'functions')));
addpath(fullfile(repository_root,'config'));
cfg = reproduction_config();

% Resolve result paths and overwrite policy; retain only Run 03 configuration.
options = struct('workflow',workflow,'saved_results_source',saved_results_source, ...
    'experiment_database_source',experiment_database_source, ...
    'save_computed_results',save_computed_results, ...
    'overwrite_generated_results',overwrite_generated_results,'figure_formats',figure_formats);
context = earth2rpo.prepare_run(repository_root,options);
scientific_cfg = struct('meta',cfg.meta,'earth_to_rpo',cfg.earth_to_rpo);

% Load the tracked orbit seeds and example indices only for recomputation.
if context.compute
    loaded = load(fullfile(repository_root,'data/inputs/run_03_reference_inputs.mat'),'inputs');
    inputs = loaded.inputs;
end
fprintf('Run 03 setup complete.\n');

%% Earth-access maps and example connections

% Recompute the backward access search or load the selected saved result.
if context.compute
    result = earth2rpo.compute_access(cfg.earth_to_rpo,inputs);

    % Attach scientific configuration and input provenance before saving.
    result = rpo.finalize_result(result,scientific_cfg,'run_03 / Earth access', ...
        'data/inputs/run_03_reference_inputs.mat','generated');
    if save_computed_results
        rpo.save_result(result,fullfile(context.result_directory,'earth_access.mat'));
    end
else
    result = rpo.load_result(fullfile(context.result_directory,'earth_access.mat'));
end

% Export the access maps/examples and report retained boundary-hit counts.
rpo2ds.export_figures(earth2rpo.plot_access(result,context.figure_visibility), ...
    context.output_directory,figure_formats);
fprintf('Earth-access counts at 100000/10000/1000 km: %d / %d / %d\n',result.metrics.candidate_counts);

%% LEO-to-RPO initial guess and optimized transfer

% Build and optimize a direct transfer, or load the selected saved solution.
if context.compute
    access = earth2rpo.load_access(repository_root,experiment_database_source);
    result = earth2rpo.compute_validation(cfg.earth_to_rpo,access.data.maps.h10000,inputs);

    % Attach the selected access-table source and retain the result if requested.
    result = rpo.finalize_result(result,scientific_cfg,'run_03 / direct LEO transfer', ...
        experiment_database_source+' Earth-access tables','generated');
    if save_computed_results
        rpo.save_result(result,fullfile(context.result_directory,'earth_to_rpo_validation.mat'));
    end
else
    result = rpo.load_result(fullfile(context.result_directory,'earth_to_rpo_validation.mat'));
end

% Export patch/node/transfer views and report physical maneuver metrics.
rpo2ds.export_figures(earth2rpo.plot_validation(result,context.figure_visibility), ...
    context.output_directory,figure_formats);
s = result.data.solution;
fprintf(['LEO altitude %.0f km; total DV %.6f m/s (departure %.6f, insertion %.6f); ', ...
    'TOF %.6f d; minimum Moon altitude %.6f km.\n'],s.departure_altitude_km, ...
    s.dv_total_mps,s.departure_summary.impulse_magnitude_mps, ...
    s.arrival_summary.impulse_magnitude_mps,s.tof_days,s.min_moon_altitude_km);

%% Departure-altitude continuation

% Solve the ascending altitude chain, or load the selected saved campaign.
if context.compute
    access = earth2rpo.load_access(repository_root,experiment_database_source);
    result = earth2rpo.compute_continuation(cfg.earth_to_rpo,access.data.maps.h10000,inputs);

    % Attach the selected access-table source and retain the campaign if requested.
    result = rpo.finalize_result(result,scientific_cfg,'run_03 / altitude continuation', ...
        experiment_database_source+' Earth-access tables','generated');
    if save_computed_results
        rpo.save_result(result,fullfile(context.result_directory,'altitude_continuation.mat'));
    end
else
    result = rpo.load_result(fullfile(context.result_directory,'altitude_continuation.mat'));
end

% Export the maneuver breakdown and display the stored altitude-sweep metrics.
rpo2ds.export_figures(earth2rpo.plot_continuation(result,context.figure_visibility), ...
    context.output_directory,figure_formats);
disp(result.data.results_table(:,{'altitude_km','total_dv_mps','departure_dv_mps', ...
    'injection_dv_mps','tof_days','min_moon_altitude_km'}));
