clc; close all; clear;
addpath(genpath('functions'));

%% Scientific purpose and instructions
% Reproduce the Comet-Interceptor-like study: warning law, heliocentric
% accessibility, RPO escape calibration and Monte Carlo comparisons
% (13 figures and probability tables).
%
% Run from the repository root. Run Setup before individual scientific sections.
% The default replays four saved results from Zenodo using base MATLAB.
% See the README for extraction instructions, or select
% saved_results_source = "local_run" to replay locally generated results.
%
% workflow = "run_experiment" recomputes with config/reproduction_config.m.
% Fresh accessibility/Monte Carlo use heliocentric_database_mode; "rebuild"
% creates databases in the accessibility section, with optional parallel workers.
% Serial RPO calibration requires Optimization Toolbox and independently
% selects its RPO-Moon table through experiment_database_source.
% Recomputation defaults to generated databases; see docs/DATA.md for
% selecting and installing reference inputs.
%
% Figures and CSV tables: outputs/{reference,generated}/06_comet_case_study/.
% Retained MAT results: results/generated/comet_case_study/.
% See README.md for shared controls and docs/REPRODUCIBILITY.md for
% capability limits, failed-scenario handling and regeneration differences.

%% User settings
workflow = "use_saved_results";          % "use_saved_results" or "run_experiment"
saved_results_source = "zenodo";            % "zenodo" or "local_run"
experiment_database_source = "generated";   % RPO-Moon: "zenodo" or "generated"
heliocentric_database_mode = "generated";   % "zenodo", "generated", or "rebuild"
save_computed_results = true;
overwrite_generated_results = false;    % Allow replacement of generated MATs/databases
figure_formats = ["png"];                % "png", "pdf", or both

%% Setup and configuration

% Resolve frozen settings and independent result/database source controls.
repository_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repository_root,'functions')));
addpath(fullfile(repository_root,'config'));
cfg = reproduction_config();
scientific_cfg = cfg.comet_case_study;
options = struct('workflow',workflow,'saved_results_source',saved_results_source, ...
    'experiment_database_source',experiment_database_source, ...
    'heliocentric_database_mode',heliocentric_database_mode, ...
    'save_computed_results',save_computed_results,'overwrite_generated_results',overwrite_generated_results);
context = comet.prepare_run(repository_root,options);
fprintf('Run06 setup complete.\n');

%% Warning-time law

% Compute the piecewise-uniform warning laws or load their saved result.
result_path = fullfile(context.result_directory,'warning_law.mat');
if context.compute
    result = comet.warning_law(scientific_cfg);
    result = rpo.finalize_result(result,cfg,mfilename,"config/reproduction_config.m","generated");
    if save_computed_results, rpo.save_result(result,result_path); end
else
    result = rpo.load_result(result_path);
end

% Export threshold/yearly probabilities and the three-panel warning-law figure.
writetable(result.threshold_table,fullfile(context.output_directory,'warning_thresholds.csv'));
writetable(result.yearly_table,fullfile(context.output_directory,'warning_yearly_probabilities.csv'));
rpo2ds.export_figures(comet.plot_warning(result,context.figure_visibility),context.output_directory,figure_formats);

%% Heliocentric databases and accessibility maps

% Optionally rebuild the Lambert grids, then extract the four study horizons.
result_path = fullfile(context.result_directory,'heliocentric_accessibility.mat');
if context.compute
    if context.rebuild, comet.rebuild_databases(scientific_cfg,context.database_directory); end
    result = comet.compute_accessibility(scientific_cfg,context.database_directory);
    inputs = [result.maps.database];
    result = rpo.finalize_result(result,cfg,mfilename,string({inputs.filename}),"generated");
    if save_computed_results, rpo.save_result(result,result_path); end
else
    result = rpo.load_result(result_path);
end

% Export all eight accessibility figures with their fixed scientific numbers.
rpo2ds.export_figures(comet.plot_accessibility(result,context.figure_visibility),context.output_directory,figure_formats);

%% RPO calibration and directional boundary

% Calibrate the physical maneuver budget and continue the planar RTH boundary.
result_path = fullfile(context.result_directory,'rpo_calibration.mat');
if context.compute
    [database,database_metadata] = rpo2ds.load_database(repository_root,experiment_database_source,'phase_aware');
    seed = rpo2ds.load_seed(repository_root);
    result = comet.compute_calibration(scientific_cfg.calibration,database,seed.state_2d);

    % Attach encounter-table identity, seed hash and explicit physical units.
    result.database = database_metadata;
    result.seed = struct('filename',"reference_rpo_seed.mat",'sha256', ...
        rpo2ds.sha256_file(fullfile(repository_root,'data','inputs','reference_rpo_seed.mat')));
    result.units = struct('impulse_cost','m/s','excess_speed','km/s','tof','days', ...
        'nodes','Earth-Moon nondimensional','angle','degrees');
    result = rpo.finalize_result(result,cfg,mfilename,["RPO-Moon phase-aware database","canonical RPO seed"],"generated");
    if save_computed_results, rpo.save_result(result,result_path); end
else
    result = rpo.load_result(result_path);
end

% Retain complete trial diagnostics in CSV alongside the directional boundary.
writetable(result.calibration.history,fullfile(context.output_directory,'rpo_calibration_history.csv'));
writetable(result.sweep,fullfile(context.output_directory,'rpo_directional_boundary.csv'));
rpo2ds.export_figures(comet.plot_calibration(result,context.figure_visibility),context.output_directory,figure_formats);

%% Monte Carlo analyses

% Evaluate both fixed architecture comparisons using one catalogue and warning sample.
result_path = fullfile(context.result_directory,'monte_carlo_analyses.mat');
if context.compute
    result = comet.compute_monte_carlo(scientific_cfg,context.database_directory);
    result = rpo.finalize_result(result,cfg,mfilename,string({result.databases.filename}),"generated");
    if save_computed_results, rpo.save_result(result,result_path); end
else
    result = rpo.load_result(result_path);
end

% Export setup, unconditional success statistics and mission-capped epoch plots.
writetable(result.combined.table,fullfile(context.output_directory,'monte_carlo_l1_l2.csv'));
writetable(result.mixed.table,fullfile(context.output_directory,'monte_carlo_mixed.csv'));
writetable(result.setup_table,fullfile(context.output_directory,'monte_carlo_setup.csv'));
disp(result.combined.table); disp(result.mixed.table);
rpo2ds.export_figures(comet.plot_monte_carlo(result,context.figure_visibility),context.output_directory,figure_formats);
