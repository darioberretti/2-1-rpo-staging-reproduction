clc; close all; clear;
addpath(genpath('functions'));

%% Scientific purpose and instructions
% Reproduce the A1 family, Moon-grazing 2:1 RPO, linear stability,
% Poincare section and vertical response, producing eight figures.
%
% Run from the repository root. Select the workflow and saved-result source
% below, then run the whole script or Setup followed by a scientific section.
% Saved-result replay uses base MATLAB and the selected installed MAT results.
% See docs/DATA.md for input locations.
%
% workflow = "run_experiment" recomputes with config/reproduction_config.m.
% Parallel Computing Toolbox is optional for the A1 sweep. Fresh computation
% takes tens of minutes and also needs reference results for comparisons.
%
% Figures: outputs/{reference,generated}/01_reference_rpo/.
% Retained MAT results: results/generated/reference_rpo/.
% See README.md for shared controls and docs/REPRODUCIBILITY.md for
% scientific interpretation.

%% User settings
workflow = "use_saved_results";       % "use_saved_results" or "run_experiment"

% Select the result source when workflow = "use_saved_results".
saved_results_source = "zenodo";      % "zenodo" or "local_run"

% Control result retention when workflow = "run_experiment".
% If false, use temporary MAT files between stages without retaining results.
save_computed_results = true;
overwrite_generated_results = false;   % Required for batch replacement

figure_formats = "png";               % "png", "pdf", or "both"

%% Setup and configuration

% Resolve the workflow, result paths, and figure export settings.
repository_root = fileparts(mfilename('fullpath'));
reference_directory = fullfile(repository_root, 'results', 'reference', 'reference_rpo');
addpath(fullfile(repository_root, 'config'));

run_options = struct( ...
    'workflow', workflow, ...
    'saved_results_source', saved_results_source, ...
    'save_computed_results', save_computed_results, ...
    'overwrite_generated_results', overwrite_generated_results, ...
    'figure_formats', figure_formats);
run_context = rpo.prepare_run(repository_root, run_options);

% Load the frozen scientific settings when recomputing the experiments.
if run_context.compute
    cfg = reproduction_config();
    if cfg.reference_rpo.a1.use_parallel && ...
            ~rpo.parallel_available(cfg.reference_rpo.a1.use_parallel)
        fprintf(['Parallel Computing Toolbox is not available. The A1 sweep ', ...
            'will run serially and take longer; install the toolbox for ', ...
            'faster execution.\n\n']);
    end
end

%% A1 family and reference RPO

% Compute the A1 reconstruction or load the selected saved result.
a1_filename = fullfile(run_context.result_directory, 'a1_reference_rpo.mat');

if run_context.compute
    result = rpo.compute_a1_reference(cfg);

    % Attach provenance and compare against the reference result.
    result = rpo.finalize_result(result, cfg, 'run_01_reference_rpo.m', ...
        "Frozen A1 sweep and physical continuation settings", ...
        run_context.artifact_type);
    result = rpo.add_reference_comparison('a1', result, fullfile( ...
        reference_directory, 'a1_reference_rpo.mat'));

    rpo.save_result(result, a1_filename);
else
    result = rpo.load_result(a1_filename);
end

% Report the reconstructed branches and reference orbit.
fprintf('A1/RPO source: %s\n', a1_filename);
for branch_name = ["k1", "k2", "k3"]
    branch = result.metrics.branch_ranges.(branch_name);
    fprintf('%s: x0 [%.2f, %.2f], C [%.6f, %.6f], T/2 [%.4f, %.4f] TU\n', ...
        branch_name, branch.x0_range, branch.jacobi_range, ...
        branch.half_period_range_tu);
end

metrics = result.metrics.rpo;
fprintf(['RPO state: [%.8f, %.8f, %.8f, %.8f]\n', ...
    'Period: %.8f TU = %.5f days; half-period: %.8f TU = %.5f days\n', ...
    'C: %.8f; Moon distance/altitude: %.1f / %.1f km; closure: %.3e\n\n'], ...
    metrics.state, metrics.period_tu, metrics.period_days, ...
    metrics.half_period_tu, metrics.half_period_days, ...
    metrics.jacobi_constant, metrics.minimum_moon_distance_km, ...
    metrics.minimum_moon_altitude_km, metrics.closure_residual);

% Create and export the A1 figures.
figures = rpo.plot_a1_reference(result, run_context.figure_visibility);
rpo.export_figures(figures, run_context.output_directory, ...
    run_context.figure_formats);

%% Linear stability

% Compute the stability analysis or load the selected saved result.
stability_filename = fullfile(run_context.result_directory, 'linear_stability.mat');

if run_context.compute
    a1_result = rpo.load_result(a1_filename);
    result = rpo.compute_stability(cfg, a1_result);

    % Attach provenance and compare against the reference result.
    result = rpo.finalize_result(result, cfg, 'run_01_reference_rpo.m', ...
        "a1_reference_rpo.mat", run_context.artifact_type);
    result = rpo.add_reference_comparison('stability', result, fullfile( ...
        reference_directory, 'linear_stability.mat'));

    rpo.save_result(result, stability_filename);
else
    result = rpo.load_result(stability_filename);
end

% Report the planar and spatial stability classifications.
fprintf('Linear-stability source: %s\n', stability_filename);

metrics = result.metrics.reference_rpo;
stability_labels = {'linearly unstable', 'linearly stable'};

fprintf('Reference RPO planar rho: %.12f (%s)\n', ...
    metrics.planar_spectral_radius, stability_labels{1 + metrics.planar_stable});
fprintf('Reference RPO spatial rho: %.12f (%s)\n\n', ...
    metrics.spatial_spectral_radius, stability_labels{1 + metrics.spatial_stable});

% Create and export the stability figures.
figures = rpo.plot_stability(result, run_context.figure_visibility);
rpo.export_figures(figures, run_context.output_directory, ...
    run_context.figure_formats);

%% Poincare section

% Compute the Poincare section or load the selected saved result.
poincare_filename = fullfile(run_context.result_directory, 'poincare_section.mat');

if run_context.compute
    a1_result = rpo.load_result(a1_filename);
    result = rpo.compute_poincare(cfg, a1_result);

    % Attach provenance and compare against the reference result.
    result = rpo.finalize_result(result, cfg, 'run_01_reference_rpo.m', ...
        "a1_reference_rpo.mat", run_context.artifact_type);
    result = rpo.add_reference_comparison('poincare', result, fullfile( ...
        reference_directory, 'poincare_section.mat'));

    rpo.save_result(result, poincare_filename);
else
    result = rpo.load_result(poincare_filename);
end

% Report the sampling settings and retained crossings.
fprintf('Poincare source: %s\n', poincare_filename);
fprintf(['Seed %d; %d requested initial conditions; %d valid chains; ', ...
    '%d crossings retained.\n\n'], ...
    result.configuration.reference_rpo.poincare.random_seed, ...
    result.metrics.requested_initial_condition_count, ...
    result.metrics.valid_chain_count, numel(result.data.crossing_x));

% Create and export the Poincare figures.
figures = rpo.plot_poincare(result, run_context.figure_visibility);
rpo.export_figures(figures, run_context.output_directory, ...
    run_context.figure_formats);

%% Vertical instability

% Compute the vertical response or load the selected saved result.
vertical_filename = fullfile(run_context.result_directory, 'vertical_instability.mat');

if run_context.compute
    a1_result = rpo.load_result(a1_filename);
    result = rpo.compute_vertical(cfg, a1_result);

    % Attach provenance and compare against the reference result.
    result = rpo.finalize_result(result, cfg, 'run_01_reference_rpo.m', ...
        "a1_reference_rpo.mat", run_context.artifact_type);
    result = rpo.add_reference_comparison('vertical', result, fullfile( ...
        reference_directory, 'vertical_instability.mat'));

    rpo.save_result(result, vertical_filename);
else
    result = rpo.load_result(vertical_filename);
end

% Report the peak vertical response across impulse magnitudes and phases.
fprintf('Vertical-response source: %s\n', vertical_filename);
for magnitude = unique(abs(result.data.magnitude.impulses_vu), 'stable')'
    rows = abs(result.data.magnitude.impulses_vu) == magnitude;
    fprintf('|vz0| = %.0e VU: maximum |z| = %.6f DU\n', magnitude, ...
        max(result.metrics.maximum_absolute_z_by_impulse_du(rows)));
end

phase_peaks = result.metrics.maximum_absolute_z_by_phase_du;
fprintf('Phase sweep at vz0 = %.0e VU: %d phases, peak |z| range [%.6f, %.6f] DU\n', ...
    result.data.phase.impulse_vu, numel(result.data.phase.fractions), ...
    min(phase_peaks), max(phase_peaks));

% Create and export the vertical-response figures.
figures = rpo.plot_vertical(result, run_context.figure_visibility);
rpo.export_figures(figures, run_context.output_directory, ...
    run_context.figure_formats);
