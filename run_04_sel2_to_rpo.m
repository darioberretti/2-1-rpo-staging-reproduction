clc; close all; clear;
addpath(genpath('functions'));

%% Scientific purpose and instructions
% Reproduce Sun-Earth L2 halo access to the RPO: reference halo, unstable
% directions, database proxy, transfer and Sun-phase continuation (nine figures).
%
% Run from the repository root. Run Setup before individual scientific sections.
% The default replays five saved results from Zenodo using base MATLAB.
% See the README for extraction instructions, or select
% saved_results_source = "local_run" to replay locally generated results.
%
% workflow = "run_experiment" recomputes with config/reproduction_config.m.
% Transfer solves require Optimization Toolbox; full continuation takes hours.
% Only the proxy section uses experiment_database_source. Transfer solves
% use tracked unoptimized shooting-node problems.
%
% Figures: outputs/{reference,generated}/04_sel2_to_rpo/.
% Retained MAT results: results/generated/sel2_to_rpo/.
% Continuation saves intermediate progress; only completed results can replay.
% See README.md for shared controls and docs/REPRODUCIBILITY.md for
% proxy interpretation and fresh-run differences.

%% User settings
workflow = "use_saved_results";        % "use_saved_results" or "run_experiment"
saved_results_source = "zenodo";        % "zenodo" or "local_run"
experiment_database_source = "zenodo"; % "zenodo" or "generated" proxy databases
save_computed_results = true;
overwrite_generated_results = false;
figure_formats = ["png"];              % "png", "pdf", or both

%% Setup and configuration

% Resolve the study configuration and saved/generated workflow paths.
repository_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repository_root,'functions')));
addpath(fullfile(repository_root,'config'));
cfg = reproduction_config();
scientific_cfg = struct('meta',cfg.meta,'sel2_to_rpo',cfg.sel2_to_rpo);
options = struct('workflow',workflow,'saved_results_source',saved_results_source, ...
    'experiment_database_source',experiment_database_source, ...
    'save_computed_results',save_computed_results,'overwrite_generated_results',overwrite_generated_results);
context = sel2rpo.prepare_run(repository_root,options);

% Load frozen unoptimized problems and branch orientation only for fresh runs.
if context.compute
    loaded = load(fullfile(repository_root,'data','inputs','run_04_reference_inputs.mat'),'inputs');
    inputs = loaded.inputs;
end
fprintf('Run 04 setup complete.\n');

%% Reference halo: analytical seed, correction and amplitude continuation
if context.compute
    % Correct and amplitude-continue a Richardson seed, then check orbit closure.
    result = sel2rpo.compute_halo(scientific_cfg);
    if save_computed_results
        rpo.save_result(result,fullfile(context.result_directory,'sel2_reference_halo.mat'));
    end
else
    result = rpo.load_result(fullfile(context.result_directory,'sel2_reference_halo.mat'));
end

% Export the phase-colored halo and report achieved physical amplitudes.
rpo2ds.export_figures(sel2rpo.plot_halo(result,context.figure_visibility),context.output_directory,figure_formats);
fprintf('Halo: period %.9f days; Ay %.3f km; Az %.3f km; closure %.3e.\n', ...
    result.metrics.period_days,result.metrics.Ay_achieved_km,result.metrics.Az_achieved_km,result.metrics.closure_inf_nd);

%% Both-branch unstable-direction illustration (x-y and x-z projections)
if context.compute
    % Recompute the halo and propagate both signs of velocity-only perturbations.
    result = sel2rpo.compute_manifolds(scientific_cfg,inputs.manifold_reference_direction);
    if save_computed_results
        rpo.save_result(result,fullfile(context.result_directory,'sel2_unstable_directions.mat'));
    end
else
    result = rpo.load_result(fullfile(context.result_directory,'sel2_unstable_directions.mat'));
end

% Export SE projections and report the dominant unstable multiplier.
rpo2ds.export_figures(sel2rpo.plot_manifolds(result,context.figure_visibility),context.output_directory,figure_formats);
fprintf('Unstable multiplier: %.9f; halo closure %.3e.\n',result.metrics.unstable_multiplier,result.metrics.closure_inf_nd);

%% Database-level proxy Pareto front (not an optimized transfer cost)
if context.compute
    % Evaluate all duration-capped encounter pairs; sample only the display cloud.
    result = sel2rpo.compute_proxy(scientific_cfg,repository_root,experiment_database_source);
    if save_computed_results
        rpo.save_result(result,fullfile(context.result_directory,'database_proxy_pareto.mat'));
    end
else
    result = rpo.load_result(fullfile(context.result_directory,'database_proxy_pareto.mat'));
end

% Export the scalar speed-mismatch frontier and report admissible pair counts.
rpo2ds.export_figures(sel2rpo.plot_proxy(result,context.figure_visibility),context.output_directory,figure_formats);
fprintf('Proxy: %d frontier rows from %.0f pairs below the duration cap.\n', ...
    result.metrics.frontier_rows,result.metrics.admissible_pairs);

%% Reconstructed initial guess and representative optimized transfer
if context.compute
    % Reconstruct the frozen zero-phase guess, optimize it and repropagate its arcs.
    result = sel2rpo.compute_validation(scientific_cfg,inputs);
    if save_computed_results
        rpo.save_result(result,fullfile(context.result_directory,'sel2_to_rpo_validation.mat'));
    end
else
    result = rpo.load_result(fullfile(context.result_directory,'sel2_to_rpo_validation.mat'));
end

% Export the construction and solved trajectories in both rotating frames.
rpo2ds.export_figures(sel2rpo.plot_validation(result,context.figure_visibility),context.output_directory,figure_formats);
fprintf('Transfer: DV %.6f m/s; waiting %.6f days; flight %.6f days; total elapsed %.6f days.\n', ...
    result.metrics.delta_v_m_s,result.metrics.waiting_days,result.metrics.flight_days,result.metrics.elapsed_days);

%% Complete historical Sun-phase continuation and campaign merge
if context.compute
    % Optionally save intermediate campaign merges and the final completed sweep.
    result_file='';
    if save_computed_results, result_file=fullfile(context.result_directory,'sel2_to_rpo_continuation.mat'); end
    result = sel2rpo.compute_continuation(scientific_cfg,inputs,result_file);
else
    result = rpo.load_result(fullfile(context.result_directory,'sel2_to_rpo_continuation.mat'));

    % Reject saved intermediate merges so a partial sweep cannot appear complete.
    if ~result.status.completed
        error('sel2rpo:IncompleteContinuation','This saved continuation contains only completed intermediate campaigns. Let the full computation finish before replaying it.');
    end
end

% Export cost and flight/elapsed times against the merged Sun-phase grid.
rpo2ds.export_figures(sel2rpo.plot_continuation(result,context.figure_visibility),context.output_directory,figure_formats);
fprintf('Continuation: %d phases; DV range %.6f--%.6f m/s.\n', ...
    result.metrics.case_count,result.metrics.delta_v_range_m_s);
