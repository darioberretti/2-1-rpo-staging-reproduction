clc; close all; clear;
addpath(genpath('functions'));

%% Scientific purpose and instructions
% Reproduce the paper's one-year stationkeeping comparison: 1000 trials per
% periodic solution, setup/results tables and one annual-cost distribution.
%
% Run from the repository root. Run Setup before the comparison section.
% Saved-result replay uses the selected installed MAT result and base MATLAB.
% See docs/DATA.md for input locations.
%
% workflow = "run_experiment" recomputes both branches serially from tracked
% inputs using config/reproduction_config.m. Computation uses base MATLAB,
% takes tens of minutes and prints progress. Existing generated results
% require overwrite_generated_results = true; there is no confirmation prompt.
%
% Figures and CSV tables: outputs/{reference,generated}/05_stationkeeping/.
% Retained MAT result: results/generated/stationkeeping/.
% See README.md for shared controls and docs/REPRODUCIBILITY.md for
% completion, targeting-tolerance and annual-cost definitions.

%% User settings
workflow = "use_saved_results";        % "use_saved_results" or "run_experiment"
saved_results_source = "zenodo";       % "zenodo" or "local_run"
save_computed_results = true;
overwrite_generated_results = false;  % Explicit permission to replace local MAT
figure_formats = ["png"];              % "png", "pdf", or both

%% Setup and configuration

% Load the frozen settings and resolve the selected result/output paths.
repository_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repository_root,'functions')));
addpath(fullfile(repository_root,'config'));
cfg = reproduction_config();
options = struct('workflow',workflow,'saved_results_source',saved_results_source, ...
    'save_computed_results',save_computed_results, ...
    'overwrite_generated_results',overwrite_generated_results);
context = stationkeeping.prepare_run(repository_root,options);

% Read corrected branch states and their provenance only for fresh computation.
if context.compute
    loaded = load(fullfile(repository_root,'data','inputs','run_05_reference_inputs.mat'),'inputs');
    inputs = loaded.inputs;
end
fprintf('Run 05 setup complete.\n');

%% One-year comparison: both branches, setup/results tables and cost CDF

% Compute both serial campaigns or load the selected self-contained result.
if context.compute
    result = stationkeeping.compute(cfg.stationkeeping,inputs);
    if save_computed_results
        rpo.save_result(result,context.result_file);
    end
else
    result = rpo.load_result(context.result_file);
end

% Export the saved experiment settings, completion counts and cost statistics.
stationkeeping.report(result,context.output_directory);
rpo2ds.export_figures(stationkeeping.plot_comparison(result,context.figure_visibility), ...
    context.output_directory,figure_formats);
