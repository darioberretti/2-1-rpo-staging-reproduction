clc; close all; clear;
addpath(genpath('functions'));

%% Scientific purpose and instructions
% Reproduce transfer-construction diagnostics, Venus/Mars validation,
% Sun-phase flexibility and directional reachability, producing 17 figures.
%
% Run from the repository root. Run Setup before individual scientific sections.
% Saved-result replay uses five installed MAT results and base MATLAB.
% See docs/DATA.md for input locations and database source selection.
%
% workflow = "run_experiment" recomputes with config/reproduction_config.m.
% Fresh validation/flexibility require CasADi/IPOPT; reachability requires
% Optimization Toolbox. Parallel workers are optional for supported sweeps.
% Frozen flexibility/reachability starts are inputs; their optimizations run.
%
% Figures: outputs/{reference,generated}/02_rpo_to_deep_space/.
% Retained MAT results: results/generated/rpo_to_deep_space/.
% See README.md for dependencies and shared controls, and
% docs/REPRODUCIBILITY.md for fresh-run agreement limits.

%% User settings
workflow = "use_saved_results";       % "use_saved_results" or "run_experiment"

% Choose the saved result source for replay.
saved_results_source = "zenodo";      % "zenodo" or "local_run"

% Choose inputs and result retention for fresh computation.
experiment_database_source = "zenodo"; % "zenodo" or "generated"
save_computed_results = true;
overwrite_generated_results = false;

figure_formats = ["png"];             % "png", "pdf", or both

%% Setup and configuration

% Locate the repository and load the frozen scientific settings.
repository_root = fileparts(mfilename('fullpath'));
if isempty(repository_root)
    repository_root = pwd;
end
addpath(genpath(fullfile(repository_root, 'functions')));
addpath(fullfile(repository_root, 'config'));
cfg = reproduction_config();

% Resolve replay/recomputation paths, overwrite policy, and figure settings.
options = struct('workflow', workflow, ...
    'saved_results_source', saved_results_source, ...
    'experiment_database_source', experiment_database_source, ...
    'save_computed_results', save_computed_results, ...
    'overwrite_generated_results', overwrite_generated_results, ...
    'figure_formats', figure_formats);
run_context = rpo2ds.prepare_run(repository_root, options);
fprintf('Run 02 setup complete.\n');

%% Transfer-construction diagnostics

% Recompute the construction diagnostics or load the selected saved result.
if run_context.compute && strcmpi(string(experiment_database_source), 'generated')
    fprintf(['The escape-direction comparison in this section used the legacy dense ', ...
        'database for the accepted result. Generated-database results are ', ...
        'marked non-reference; differences are expected.\n']);
end
if run_context.compute
    targets = rpo2ds.load_targets(repository_root);
    [phase_database, phase_metadata] = rpo2ds.load_database(repository_root, ...
        experiment_database_source, 'phase_aware');
    if strcmpi(string(experiment_database_source), 'generated')
        % The public database builder creates the phase-aware table. Reusing
        % it here keeps one generated database sufficient for every section.
        altitude_database = phase_database;
        altitude_metadata = phase_metadata;
    else
        [altitude_database, altitude_metadata] = rpo2ds.load_database( ...
            repository_root, experiment_database_source, 'legacy_dense');
    end
    result = rpo2ds.compute_method_exploration(cfg, phase_database, ...
        altitude_database, targets);

    % Attach the frozen configuration and input provenance.
    result = rpo2ds.finalize_result(result, cfg, ...
        'run_02_rpo_to_deep_space / method exploration', ...
        run_context.artifact_type, struct('phase_aware_database', phase_metadata, ...
        'altitude_database', altitude_metadata, 'inputs', ...
        rpo2ds.input_metadata(repository_root, true)));

    % Retain the compact result when requested and compare against accepted data.
    if save_computed_results
        rpo2ds.save_result(result, fullfile(run_context.result_directory, ...
            'method_exploration.mat'), true);
    end
    rpo2ds.compare_result(result, fullfile(repository_root, 'results', 'reference', ...
        'rpo_to_deep_space', 'method_exploration.mat'), 'method exploration', ...
        comparison_file(run_context, save_computed_results, 'method_exploration'));
else
    result = rpo2ds.load_result(fullfile(run_context.result_directory, ...
        'method_exploration.mat'));
end

% Create and export this section's figures under their fixed scientific numbers.
figures = rpo2ds.plot_method_exploration(result, run_context.figure_visibility);
rpo2ds.export_figures(figures, run_context.output_directory, run_context.figure_formats);
% Report the stored Earth-relative targets used by the transfer construction.
if isfield(result, 'metrics') && isfield(result.metrics, 'targets')
    fprintf('Accepted Earth-relative target table:\n');
    for target_name = ["venus", "mars"]
        target_metrics = result.metrics.targets.(target_name);
        fprintf(['  %-6s |v_inf| %.6f km/s, TOF %.6f d, ECI [% .6f % .6f % .6f] ', ...
            'RTH [% .6f % .6f % .6f] km/s\n'], ...
            char(upper(target_name)), target_metrics.vinf_earth_magnitude_km_s, ...
            target_metrics.time_of_flight_days, target_metrics.vinf_earth_eci_km_s, ...
            target_metrics.vinf_earth_rth_km_s);
    end
end
fprintf('Diagnostics complete: database rows and accepted target summaries are in method_exploration.mat.\n');

%% Venus validation

% Solve the Venus validation chain or load the selected saved result.
if run_context.compute
    targets = rpo2ds.load_targets(repository_root);
    seed = rpo2ds.load_seed(repository_root);
    [phase_database, phase_metadata] = rpo2ds.load_database(repository_root, ...
        experiment_database_source, 'phase_aware');
    result = rpo2ds.compute_validation(cfg, phase_database, seed, 'venus', targets.venus);

    % Attach the frozen configuration and input provenance.
    result = rpo2ds.finalize_result(result, cfg, ...
        'run_02_rpo_to_deep_space / Venus validation', run_context.artifact_type, ...
        struct('database', phase_metadata, 'inputs', ...
        rpo2ds.input_metadata(repository_root, true)));

    % Retain the compact result when requested and compare against accepted data.
    if save_computed_results
        rpo2ds.save_result(result, fullfile(run_context.result_directory, ...
            'venus_validation.mat'), true);
    end
    rpo2ds.compare_result(result, fullfile(repository_root, 'results', 'reference', ...
        'rpo_to_deep_space', 'venus_validation.mat'), 'Venus validation', ...
        comparison_file(run_context, save_computed_results, 'venus_validation'));
else
    result = rpo2ds.load_result(fullfile(run_context.result_directory, ...
        'venus_validation.mat'));
end

% Create and export this section's figures under their fixed scientific numbers.
figures = rpo2ds.plot_validation(result, 'venus', run_context.figure_visibility);
rpo2ds.export_figures(figures, run_context.output_directory, run_context.figure_formats);

% Report the stored section metrics.
rpo2ds.print_summary(result, 'venus');
fprintf('Venus validation complete.\n');

%% Mars validation

% Solve the Mars validation chain or load the selected saved result.
if run_context.compute
    targets = rpo2ds.load_targets(repository_root);
    seed = rpo2ds.load_seed(repository_root);
    [phase_database, phase_metadata] = rpo2ds.load_database(repository_root, ...
        experiment_database_source, 'phase_aware');
    result = rpo2ds.compute_validation(cfg, phase_database, seed, 'mars', targets.mars);

    % Attach the frozen configuration and input provenance.
    result = rpo2ds.finalize_result(result, cfg, ...
        'run_02_rpo_to_deep_space / Mars validation', run_context.artifact_type, ...
        struct('database', phase_metadata, 'inputs', ...
        rpo2ds.input_metadata(repository_root, true)));

    % Retain the compact result when requested and compare against accepted data.
    if save_computed_results
        rpo2ds.save_result(result, fullfile(run_context.result_directory, ...
            'mars_validation.mat'), true);
    end
    rpo2ds.compare_result(result, fullfile(repository_root, 'results', 'reference', ...
        'rpo_to_deep_space', 'mars_validation.mat'), 'Mars validation', ...
        comparison_file(run_context, save_computed_results, 'mars_validation'));
else
    result = rpo2ds.load_result(fullfile(run_context.result_directory, ...
        'mars_validation.mat'));
end

% Create and export this section's figures under their fixed scientific numbers.
figures = rpo2ds.plot_validation(result, 'mars', run_context.figure_visibility);
rpo2ds.export_figures(figures, run_context.output_directory, run_context.figure_formats);

% Report the stored section metrics.
rpo2ds.print_summary(result, 'mars');
fprintf('Mars validation complete.\n');

%% Sun-phase flexibility

% Run the frozen Sun-phase campaigns or load the selected saved result.
if run_context.compute
    targets = rpo2ds.load_targets(repository_root);
    [flexibility_preprocessing, flexibility_preprocessing_metadata] = ...
        rpo2ds.load_flexibility_initial_guesses(repository_root);
    [phase_database, phase_metadata] = rpo2ds.load_database(repository_root, ...
        experiment_database_source, 'phase_aware');
    result = rpo2ds.compute_flexibility(cfg, phase_database, targets, ...
        flexibility_preprocessing);

    % Attach the frozen configuration and input provenance.
    result = rpo2ds.finalize_result(result, cfg, ...
        'run_02_rpo_to_deep_space / Sun-phase flexibility', run_context.artifact_type, ...
        struct('database', phase_metadata, ...
        'flexibility_preprocessing', flexibility_preprocessing_metadata, 'inputs', ...
        rpo2ds.input_metadata(repository_root, true)));

    % Retain the compact result when requested and compare against accepted data.
    if save_computed_results
        rpo2ds.save_result(result, fullfile(run_context.result_directory, ...
            'flexibility.mat'), true);
    end
    rpo2ds.compare_result(result, fullfile(repository_root, 'results', 'reference', ...
        'rpo_to_deep_space', 'flexibility.mat'), 'Sun-phase flexibility', ...
        comparison_file(run_context, save_computed_results, 'flexibility'));
else
    result = rpo2ds.load_result(fullfile(run_context.result_directory, ...
        'flexibility.mat'));
end

% Create and export this section's figures under their fixed scientific numbers.
figures = rpo2ds.plot_flexibility(result, run_context.figure_visibility);
rpo2ds.export_figures(figures, run_context.output_directory, run_context.figure_formats);

% Report the stored section metrics.
rpo2ds.print_summary(result, 'flexibility');
fprintf('Sun-phase flexibility complete.\n');

%% Directional reachability

% Run the historical directional sweep or load the selected saved result.
if run_context.compute && strcmpi(string(experiment_database_source), 'generated')
    fprintf(['The accepted reachability evidence used the legacy dense database. ', ...
        'Generated-database results are marked non-reference; differences are expected.\n']);
end
if run_context.compute
    [reachability_preprocessing, preprocessing_metadata] = ...
        rpo2ds.load_reachability_initial_guesses(repository_root);
    if strcmpi(string(experiment_database_source), 'generated')
        % One generated phase-aware table serves every fresh section.
        [reachability_database, reachability_metadata] = rpo2ds.load_database( ...
            repository_root, experiment_database_source, 'phase_aware');
    else
        % The accepted directional sweep used the historical dense table.
        [reachability_database, reachability_metadata] = rpo2ds.load_database( ...
            repository_root, experiment_database_source, 'legacy_dense');
    end
    result = rpo2ds.compute_reachability(cfg, reachability_database, ...
        reachability_preprocessing);

    % Attach the frozen configuration and input provenance.
    result = rpo2ds.finalize_result(result, cfg, ...
        'run_02_rpo_to_deep_space / directional reachability', run_context.artifact_type, ...
        struct('database', reachability_metadata, ...
        'reachability_preprocessing', preprocessing_metadata, 'inputs', ...
        rpo2ds.input_metadata(repository_root, false)));

    % Retain the compact result when requested and compare against accepted data.
    if save_computed_results
        rpo2ds.save_result(result, fullfile(run_context.result_directory, ...
            'reachability.mat'), true);
    end
    rpo2ds.compare_result(result, fullfile(repository_root, 'results', 'reference', ...
        'rpo_to_deep_space', 'reachability.mat'), 'directional reachability', ...
        comparison_file(run_context, save_computed_results, 'reachability'));
else
    result = rpo2ds.load_result(fullfile(run_context.result_directory, ...
        'reachability.mat'));
end

% Create and export this section's figures under their fixed scientific numbers.
figures = rpo2ds.plot_reachability(result, run_context.figure_visibility);
rpo2ds.export_figures(figures, run_context.output_directory, run_context.figure_formats);

% Report the stored section metrics.
rpo2ds.print_summary(result, 'reachability');
fprintf('Directional reachability complete.\n');

function filename = comparison_file(run_context, save_results, result_name)
%COMPARISON_FILE Select a comparison MAT destination only when saving results.

if save_results
    filename = fullfile(run_context.result_directory, [result_name '_comparison.mat']);
else
    filename = [];
end
end
