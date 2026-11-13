function context = prepare_run(repository_root,options)
%PREPARE_RUN Select saved/generated paths and check result replacement.
%
% PROTOTYPE:
%   context = sel2rpo.prepare_run(repository_root, options)
%
% DESCRIPTION:
%   Resolve workflow and result source before any scientific section.
%   Computation uses generated results; replay selects reference for
%   zenodo or generated for local_run. Existing generated MAT products
%   require desktop confirmation unless overwrite is enabled; batch
%   replacement requires the explicit overwrite flag.
%
% INPUT:
%   repository_root  Repository directory path.
%   options          Struct with workflow (use_saved_results or
%                      run_experiment), saved_results_source (zenodo or
%                      local_run), experiment_database_source (zenodo or
%                      generated), save_computed_results and
%                      overwrite_generated_results (logical flags).
%
% OUTPUT:
%   context          Struct with compute (logical), figure_visibility
%                      (on/off), result_directory and output_directory.
%                      Paths are selected here; directories are not created.

% Resolve the workflow, display mode and selected result branch.
context.compute = strcmp(validatestring(char(options.workflow), ...
    {'use_saved_results','run_experiment'}),'run_experiment');
context.figure_visibility='off';
if usejava('desktop'), context.figure_visibility='on'; end
if context.compute
    validatestring(char(options.experiment_database_source),{'zenodo','generated'});
    branch='generated';
else
    source=validatestring(char(options.saved_results_source),{'zenodo','local_run'});
    branch='reference'; if strcmp(source,'local_run'), branch='generated'; end
end
context.result_directory=fullfile(repository_root,'results',branch,'sel2_to_rpo');
context.output_directory=fullfile(repository_root,'outputs',branch,'04_sel2_to_rpo');

% Check all five products before starting an experiment that will save results.
if context.compute && options.save_computed_results && ~options.overwrite_generated_results
    names=["sel2_reference_halo.mat","sel2_unstable_directions.mat","database_proxy_pareto.mat", ...
        "sel2_to_rpo_validation.mat","sel2_to_rpo_continuation.mat"];
    existing=names(arrayfun(@(name) isfile(fullfile(context.result_directory,name)),names));
    if ~isempty(existing)
        if ~usejava('desktop')
            error('sel2rpo:ResultsExist','Generated Run 04 results exist. Set overwrite_generated_results = true for batch replacement.');
        end
        fprintf('Generated Run 04 files to replace in %s:\n',context.result_directory);
        fprintf('  %s\n',existing);
        answer=lower(strtrim(input('Proceed? [y/N]: ','s')));
        if ~ismember(answer,{'y','yes'})
            error('sel2rpo:OverwriteCancelled','Generated-result replacement cancelled.');
        end
    end
end

% Report the resolved paths for both full runs and section-only use.
fprintf('Run 04 workflow: %s; result branch: %s\n',options.workflow,branch);
fprintf('Results: %s\nFigures: %s\n\n',context.result_directory,context.output_directory);
end
