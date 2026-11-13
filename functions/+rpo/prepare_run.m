function context = prepare_run(repository_root, options)
%PREPARE_RUN Resolve result/output paths and protect generated result slots.
%
% PROTOTYPE:
%   context = rpo.prepare_run(repository_root, options)
%
% DESCRIPTION:
%   Resolve run01 replay or recomputation paths, normalize export options,
%   and create required output directories. Recomputations require the frozen
%   reference files for comparison. Existing generated results require either
%   explicit overwrite permission or interactive confirmation.
%
% INPUT:
%   repository_root  Repository root path.
%   options          Struct with workflow (use_saved_results/run_experiment),
%                      saved_results_source (zenodo/local_run), logical
%                      save_computed_results / overwrite_generated_results,
%                      and figure_formats (png/pdf/both).
%
% OUTPUT:
%   context          Struct with compute, artifact_type, result_directory,
%                      output_directory, result_filenames, figure_formats,
%                      and figure_visibility. When saving is disabled, a
%                      transient result directory is owned by the returned
%                      onCleanup handle; keep context alive while using it.

% Normalize the workflow and common figure settings.
workflow = validatestring(options.workflow, ...
    {'use_saved_results', 'run_experiment'});
context.repository_root = repository_root;
context.compute = strcmp(workflow, 'run_experiment');
context.result_filenames = ["a1_reference_rpo.mat", "linear_stability.mat", ...
    "poincare_section.mat", "vertical_instability.mat"];
context.figure_formats = local_figure_formats(options.figure_formats);
context.figure_visibility = local_figure_visibility();
context.transient_results = false;

% Select recomputation or saved-result paths and enforce the result policy.
if context.compute
    context.artifact_type = "generated";
    context.result_directory = fullfile(repository_root, 'results', ...
        'generated', 'reference_rpo');
    context.output_directory = fullfile(repository_root, 'outputs', ...
        'generated', '01_reference_rpo');
    local_require_reference_results(repository_root, context.result_filenames);

    % Keep intermediate MAT files temporary when reusable results are not requested.
    if ~options.save_computed_results
        context.result_directory = tempname;
        context.transient_results = true;
        transient_directory = context.result_directory;
        context.transient_cleanup = onCleanup( ...
            @() local_remove_transient(transient_directory));
    else
        local_confirm_overwrite(context.result_directory, ...
            context.result_filenames, options.overwrite_generated_results);
    end
    fprintf('Run 01 will recompute the frozen experiment.\n');
    fprintf('Completed results: %s\n', context.result_directory);
else
    source = validatestring(options.saved_results_source, ...
        {'zenodo', 'local_run'});
    switch source
        case 'zenodo'
            context.artifact_type = "reference";
            context.result_directory = fullfile(repository_root, 'results', ...
                'reference', 'reference_rpo');
            context.output_directory = fullfile(repository_root, 'outputs', ...
                'reference', '01_reference_rpo');
        case 'local_run'
            context.artifact_type = "generated";
            context.result_directory = fullfile(repository_root, 'results', ...
                'generated', 'reference_rpo');
            context.output_directory = fullfile(repository_root, 'outputs', ...
                'generated', '01_reference_rpo');
    end
    fprintf('Run 01 will recreate outputs from %s saved results.\n', source);
    fprintf('Result source: %s\n', context.result_directory);
end

% Create only the directories that this workflow writes to.
if context.compute && ~isfolder(context.result_directory)
    mkdir(context.result_directory);
end
if ~isfolder(context.output_directory)
    mkdir(context.output_directory);
end
fprintf('Figure output: %s\n\n', context.output_directory);
end

function local_remove_transient(directory)
%LOCAL_REMOVE_TRANSIENT Remove the temporary result directory on context cleanup.
%
% The directory must exist and begin with the system temporary-directory path.

if isfolder(directory) && startsWith(directory, tempdir, 'IgnoreCase', true)
    rmdir(directory, 's');
end
end

function local_require_reference_results(repository_root, filenames)
%LOCAL_REQUIRE_REFERENCE_RESULTS Require all frozen MAT files used for comparison.
%
% repository_root locates results/reference/reference_rpo; filenames is
% the requested result-name array. Missing files produce an installation hint.

directory = fullfile(repository_root, 'results', 'reference', 'reference_rpo');
missing = filenames(~arrayfun(@(name) isfile(fullfile(directory, name)), filenames));
if isempty(missing)
    return
end
error('rpo:ReferenceResultsMissing', [ ...
    'Frozen recomputation requires these reference result files for comparison:\n', ...
    '  %s\n\nExtract accepted_results.zip into the repository root (see README.md).'], ...
    strjoin(missing, '\n  '));
end

function formats = local_figure_formats(requested)
%LOCAL_FIGURE_FORMATS Normalize png/pdf/both to a unique png/pdf string array.

formats = lower(string(requested));
if isscalar(formats) && formats == "both"
    formats = ["png", "pdf"];
end
if isempty(formats) || any(~ismember(formats, ["png", "pdf"]))
    error('rpo:FigureFormats', ...
        'figure_formats must contain "png", "pdf", or "both".');
end
formats = unique(formats, 'stable');
end

function visibility = local_figure_visibility()
%LOCAL_FIGURE_VISIBILITY Show figures in desktop MATLAB and hide batch figures.

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
%   Return when no requested files exist or overwrite is authorized.
%   Otherwise require interactive desktop confirmation; batch mode fails
%   with an instruction to enable replacement explicitly.
%
% INPUT:
%   directory          Generated result directory.
%   filenames          Result filenames to check within directory.
%   overwrite_allowed  Logical permitting replacement without prompting.
%
% OUTPUT:
%   None. Raises an error if overwrite is declined or unavailable.

existing = filenames(arrayfun(@(name) isfile(fullfile(directory, name)), filenames));
if isempty(existing)
    return
end
if overwrite_allowed
    return
end
if ~usejava('desktop')
    error('rpo:GeneratedResultsExist', [ ...
        'Generated results already exist in:\n%s\n', ...
        'Set overwrite_generated_results = true for batch replacement.'], ...
        directory);
end
fprintf('The following generated results will be replaced in %s:\n', directory);
fprintf('  %s\n', existing);
answer = lower(strtrim(input('Proceed? [y/N]: ', 's')));
if ~ismember(answer, {'y', 'yes'})
    error('rpo:OverwriteCancelled', 'Generated-result replacement was cancelled.');
end
end
