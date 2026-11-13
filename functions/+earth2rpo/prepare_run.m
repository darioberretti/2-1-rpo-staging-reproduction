function context = prepare_run(repository_root, options)
%PREPARE_RUN Resolve Run 03 paths and protect existing generated results.
%
% PROTOTYPE:
%   context = earth2rpo.prepare_run(repository_root, options)
%
% DESCRIPTION:
%   Select accepted/local replay or fresh computation, set desktop/batch
%   figure visibility, and require overwrite authorization for retained
%   generated files. Scientific data and solver dependencies are loaded later.
%
% INPUT:
%   repository_root  Repository root path.
%   options          Struct with workflow (use_saved_results/run_experiment),
%                      saved_results_source (zenodo/local_run),
%                      experiment_database_source (zenodo/generated), and
%                      logical save_computed_results/overwrite_generated_results.
%
% OUTPUT:
%   context          Struct with compute, figure_visibility, result_directory,
%                      and output_directory. Paths select reference/generated
%                      earth_to_rpo results and 03_earth_to_rpo figures.

% Select replay or recomputation and desktop/batch figure visibility.
context.compute = strcmp(validatestring(char(options.workflow), ...
    {'use_saved_results', 'run_experiment'}), 'run_experiment');
context.figure_visibility = 'off';
if usejava('desktop'), context.figure_visibility = 'on'; end
if context.compute
    validatestring(char(options.experiment_database_source), {'zenodo','generated'});
    branch = 'generated';
else
    source = validatestring(char(options.saved_results_source), {'zenodo','local_run'});
    branch = 'reference';
    if strcmp(source, 'local_run'), branch = 'generated'; end
end

% Resolve the result/figure destinations for the selected artifact branch.
context.result_directory = fullfile(repository_root, 'results', branch, 'earth_to_rpo');
context.output_directory = fullfile(repository_root, 'outputs', branch, '03_earth_to_rpo');

% Protect existing retained generated files before the expensive computations.
if context.compute && options.save_computed_results && ~options.overwrite_generated_results
    names = ["earth_access.mat", "earth_to_rpo_validation.mat", "altitude_continuation.mat"];
    existing = names(arrayfun(@(name) isfile(fullfile(context.result_directory, name)), names));
    if ~isempty(existing)
        if ~usejava('desktop')
            error('earth2rpo:ResultsExist', ...
                'Generated Run 03 results exist. Set overwrite_generated_results = true for batch replacement.');
        end
        fprintf('Generated Run 03 files to replace in %s:\n', context.result_directory);
        fprintf('  %s\n', existing);
        answer = lower(strtrim(input('Proceed? [y/N]: ', 's')));
        if ~ismember(answer, {'y','yes'})
            error('earth2rpo:OverwriteCancelled', 'Generated-result replacement cancelled.');
        end
    end
end
fprintf('Run 03 workflow: %s; result branch: %s\n', options.workflow, branch);
if context.compute, fprintf('Earth-access input source: %s\n', options.experiment_database_source); end
fprintf('Results: %s\nFigures: %s\n\n', context.result_directory, context.output_directory);
end
