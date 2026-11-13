function context = prepare_run(repository_root,options)
%PREPARE_RUN Select the stationkeeping result slot and enforce explicit replacement.
%
% PROTOTYPE:
%   context = stationkeeping.prepare_run(repository_root, options)
%
% DESCRIPTION:
%   Fresh runs use generated results. Saved replay chooses
%   reference for zenodo or generated for local_run. Check the single
%   result slot before computation; replacing it requires the explicit
%   overwrite flag in both desktop and batch execution.
%
% INPUT:
%   repository_root  Repository path.
%   options          Struct with workflow (run_experiment or
%                      use_saved_results), saved_results_source (zenodo
%                      or local_run), save_computed_results and
%                      overwrite_generated_results (logical flags).
%                      Workflow/source settings use scalar strings.
%
% OUTPUT:
%   context          Struct with compute (logical), result_file
%                      (stationkeeping_comparison.mat), output_directory
%                      and figure_visibility (on/off). This selects paths;
%                      saving/exporting creates directories later.

% Resolve the explicitly selected workflow and reference/generated branch.
context.compute = options.workflow == "run_experiment";
if context.compute
    branch = 'generated';
elseif options.workflow == "use_saved_results"
    switch options.saved_results_source
        case "zenodo", branch = 'reference';
        case "local_run", branch = 'generated';
        otherwise, error('stationkeeping:Source','Use saved_results_source = "zenodo" or "local_run".');
    end
else
    error('stationkeeping:Workflow','Use workflow = "use_saved_results" or "run_experiment".');
end

% Keep one MAT product and one numbered figure directory for each branch.
context.result_file = fullfile(repository_root,'results',branch, ...
    'stationkeeping','stationkeeping_comparison.mat');
context.output_directory = fullfile(repository_root,'outputs',branch,'05_stationkeeping');
context.figure_visibility = 'off';
if usejava('desktop'), context.figure_visibility = 'on'; end
% Catch an occupied output slot before spending time on the full campaign.
if context.compute && options.save_computed_results && ...
        ~options.overwrite_generated_results && isfile(context.result_file)
    error('stationkeeping:Overwrite', ...
        'Generated result already exists:\n%s\nSet overwrite_generated_results = true to replace it.', ...
        context.result_file);
end
fprintf('Run 05 workflow: %s; result branch: %s\nResult: %s\nFigures: %s\n\n', ...
    options.workflow,branch,context.result_file,context.output_directory);
end
