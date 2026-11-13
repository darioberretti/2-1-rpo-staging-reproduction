function context = prepare_run(repository_root, options)
%PREPARE_RUN Resolve Run 02 artifact paths and protect generated results.
%
% PROTOTYPE:
%   context = rpo2ds.prepare_run(repository_root, options)
%
% DESCRIPTION:
%   Select accepted/local replay or fresh computation, normalize export
%   formats, and create the writable output folders. This does not load
%   scientific inputs or require a numerical solver.
%
% INPUT:
%   repository_root  Repository root path.
%   options          Struct with workflow (use_saved_results/run_experiment),
%                      saved_results_source (zenodo/local_run),
%                      experiment_database_source (zenodo/generated), logical
%                      save_computed_results / overwrite_generated_results,
%                      and figure_formats (png/pdf/both).
%
% OUTPUT:
%   context          Struct with compute, artifact_type, result/output
%                      directories, result_filenames, normalized figure_formats,
%                      and figure_visibility. Saving computed MAT results is
%                      controlled by the runner rather than temporary files.

% Normalize the workflow and figure export settings.
context.repository_root = repository_root;
workflow = validatestring(char(options.workflow), {'use_saved_results', 'run_experiment'});
context.compute = strcmp(workflow, 'run_experiment');
context.result_filenames = [ ...
    "method_exploration.mat", "venus_validation.mat", "mars_validation.mat", ...
    "flexibility.mat", "reachability.mat"];
context.figure_formats = local_figure_formats(options.figure_formats);
context.figure_visibility = local_figure_visibility();

% Select the result branch and protect retained generated artifacts.
if context.compute
    validatestring(char(options.experiment_database_source), {'zenodo', 'generated'});
    context.artifact_type = "generated";
    context.result_directory = fullfile(repository_root, 'results', ...
        'generated', 'rpo_to_deep_space');
    context.output_directory = fullfile(repository_root, 'outputs', ...
        'generated', '02_rpo_to_deep_space');
    if options.save_computed_results
        local_confirm_overwrite(context.result_directory, ...
            context.result_filenames, options.overwrite_generated_results);
    end
    fprintf('Run 02 will recompute the frozen RPO-to-deep-space experiments.\n');
    fprintf('Database source: %s\n', options.experiment_database_source);
else
    source = validatestring(char(options.saved_results_source), ...
        {'zenodo', 'local_run'});
    if strcmp(source, 'zenodo')
        context.artifact_type = "reference";
        branch = 'reference';
    else
        context.artifact_type = "generated";
        branch = fullfile('generated');
    end
    context.result_directory = fullfile(repository_root, 'results', branch, ...
        'rpo_to_deep_space');
    context.output_directory = fullfile(repository_root, 'outputs', branch, ...
        '02_rpo_to_deep_space');
    fprintf('Run 02 will recreate outputs from %s saved results.\n', source);
end

% Create only the folders written by the selected workflow.
if ~isfolder(context.result_directory) && context.compute && options.save_computed_results
    mkdir(context.result_directory);
end
if ~isfolder(context.output_directory)
    mkdir(context.output_directory);
end
fprintf('Result source/destination: %s\n', context.result_directory);
fprintf('Figure output: %s\n\n', context.output_directory);
end

function formats = local_figure_formats(requested)
%LOCAL_FIGURE_FORMATS Normalize png/pdf/both to a unique png/pdf string array.

formats = lower(string(requested));
if isscalar(formats) && formats == "both"
    formats = ["png", "pdf"];
end
if isempty(formats) || any(~ismember(formats, ["png", "pdf"]))
    error('rpo2ds:FigureFormats', ...
        'figure_formats must contain "png", "pdf", or "both".');
end
formats = unique(formats, 'stable');
end

function visibility = local_figure_visibility()
%LOCAL_FIGURE_VISIBILITY Show desktop figures and hide figures for batch execution.

if usejava('desktop')
    visibility = 'on';
else
    visibility = 'off';
end
end

function local_confirm_overwrite(directory, filenames, overwrite_allowed)
%LOCAL_CONFIRM_OVERWRITE Protect existing generated result files.
%
% PROTOTYPE:
%   local_confirm_overwrite(directory, filenames, overwrite_allowed)
%
% DESCRIPTION:
%   Return when no requested file exists or overwrite is explicitly enabled.
%   Otherwise require desktop confirmation; batch execution fails with an
%   instruction to authorize replacement in the runner settings.
%
% INPUT:
%   directory          Generated-result directory.
%   filenames          Result filenames checked inside that directory.
%   overwrite_allowed  Logical allowing replacement without an interactive prompt.
%
% OUTPUT:
%   None. Raises an error if overwrite is unavailable or declined.

existing = filenames(arrayfun(@(name) isfile(fullfile(directory, name)), filenames));
if isempty(existing) || overwrite_allowed
    return
end
if ~usejava('desktop')
    error('rpo2ds:GeneratedResultsExist', [ ...
        'Generated Run 02 results already exist in:\n%s\n', ...
        'Set overwrite_generated_results = true for batch replacement.'], ...
        directory);
end
fprintf('The following generated Run 02 results will be replaced in %s:\n', directory);
fprintf('  %s\n', existing);
answer = lower(strtrim(input('Proceed? [y/N]: ', 's')));
if ~ismember(answer, {'y', 'yes'})
    error('rpo2ds:OverwriteCancelled', ...
        'Generated-result replacement was cancelled.');
end
end
