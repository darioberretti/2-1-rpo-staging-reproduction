function context = prepare_run(repository_root,options)
%PREPARE_RUN Resolve result/database sources and check generated-file replacement.
%
% PROTOTYPE:
%   context = comet.prepare_run(repository_root, options)
%
% DESCRIPTION:
%   Keep saved-result replay separate from fresh-computation data
%   controls. Replay reads reference (zenodo) or generated (local_run)
%   results without consulting database settings. Fresh computation
%   uses generated results and the selected heliocentric data branch;
%   rebuild selects generated databases and enables their replacement.
%
% INPUT:
%   repository_root  Repository path.
%   options          Struct with workflow (use_saved_results or
%                      run_experiment), saved_results_source (zenodo or
%                      local_run), experiment_database_source (zenodo or
%                      generated RPO-Moon tables), heliocentric_database_mode
%                      (zenodo, generated or rebuild), save_computed_results
%                      and overwrite_generated_results (logical flags).
%
% OUTPUT:
%   context          Struct with compute/rebuild flags, figure_visibility
%                      (on/off), result/output/database_directory paths
%                      and artifact_type (reference/generated). Creates
%                      the output folder. Existing selected MAT targets
%                      require desktop confirmation unless overwrite is
%                      enabled; batch replacement requires the flag.

% Resolve replay or computation before interpreting the database controls.
workflow = validatestring(char(options.workflow),{'use_saved_results','run_experiment'});
context.compute = strcmp(workflow,'run_experiment');
context.figure_visibility = 'off';
if usejava('desktop'), context.figure_visibility = 'on'; end
if context.compute
    branch = 'generated';
    validatestring(char(options.experiment_database_source),{'zenodo','generated'});
    mode = validatestring(char(options.heliocentric_database_mode),{'zenodo','generated','rebuild'});
else
    source = validatestring(char(options.saved_results_source),{'zenodo','local_run'});
    branch = 'generated';
    if strcmp(source,'zenodo'), branch = 'reference'; end
    mode = 'zenodo'; % Replay never consults database controls.
end

% Resolve independent result, figure and heliocentric database paths.
database_branch = 'generated';
if strcmp(mode,'zenodo'), database_branch = 'reference'; end
context.result_directory = fullfile(repository_root,'results',branch,'comet_case_study');
context.output_directory = fullfile(repository_root,'outputs',branch,'06_comet_case_study');
context.database_directory = fullfile(repository_root,'data',database_branch,'comet_lambert');
context.rebuild = context.compute && strcmp(mode,'rebuild');
context.artifact_type = string(branch);

% One early confirmation protects all files the selected computation may replace.
targets = strings(0,1);
if context.compute && options.save_computed_results
    names = ["warning_law.mat","heliocentric_accessibility.mat", ...
        "rpo_calibration.mat","monte_carlo_analyses.mat"];
    targets = fullfile(context.result_directory,names(:));
end
if context.rebuild
    targets = [targets;fullfile(context.database_directory,comet.database_file())];
end

% Check every selected result/rebuild target before starting expensive work.
existing = targets(arrayfun(@isfile,targets));
if ~isempty(existing) && ~options.overwrite_generated_results
    if ~usejava('desktop')
        error('comet:Overwrite','Generated files already exist. Set overwrite_generated_results=true for batch replacement.\n%s',join(existing,newline));
    end
    fprintf('Run06 may replace these generated files:\n');
    fprintf('  %s\n',existing);
    answer = lower(strtrim(input('Proceed? [y/N]: ','s')));
    if ~ismember(answer,{'y','yes'})
        error('comet:Cancelled','Generated-file replacement cancelled.');
    end
end

% Prepare the figure/table destination and report the selected sources.
if ~isfolder(context.output_directory), mkdir(context.output_directory); end
fprintf('Run06 workflow: %s; result branch: %s\n',workflow,branch);
if context.compute
    fprintf('Heliocentric databases: %s; RPO-Moon source: %s\n',mode,options.experiment_database_source);
end
fprintf('Results: %s\nFigures: %s\n',context.result_directory,context.output_directory);
end
