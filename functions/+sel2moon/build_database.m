function result_explorer = build_database(settings, execution, database_file, checkpoint_file)
%BUILD_DATABASE Build the SEL2 lunar-entry table with resumable ordered batches.
%
% PROTOTYPE:
%   result_explorer = sel2moon.build_database(settings, execution, ...
%       database_file, checkpoint_file)
%
% DESCRIPTION:
%   Reconstruct the halo and its unstable direction, then propagate all
%   plus-branch seeds. Keep encounters in seed/event order, without dense
%   histories. Resume compatible scientific progress and replace a completed
%   database only after its new version has been saved successfully.
%
% INPUT:
%   settings         Scientific profile defined in build_sel2_to_moon_database.m.
%   execution        Struct: use_parallel [logical], batch_size [positive
%                    integer], overwrite [logical]. These controls do not
%                    change the scientific grid or encounter ordering.
%   database_file    Generated MAT destination; reference files are not used.
%   checkpoint_file  Incomplete-progress MAT path, separate from the output.
%
% OUTPUT:
%   result_explorer  Struct with SEL2MoonDatabase [H x 22 table], settings,
%                    orbit, unstable-direction information, diagnostics and
%                    provenance. H is the number of Moon-SOI entries. State,
%                    time and velocity units are described in compute_seed.

% Protect completed data before reconstructing the halo or starting workers.
if isfile(database_file) && ~execution.overwrite
    if ~usejava('desktop')
        error('sel2moon:DatabaseExists', ...
            ['Generated database exists: %s\n' ...
             'Set overwrite_generated_database = true to replace it.'], database_file);
    end
    fprintf('The generated database will be replaced after successful completion:\n%s\n', ...
        database_file);
    answer = lower(strtrim(input('Proceed? [y/N]: ', 's')));
    if ~ismember(answer, {'y', 'yes'})
        error('sel2moon:Cancelled', 'Database replacement cancelled.');
    end
end
if execution.batch_size < 1 || execution.batch_size ~= fix(execution.batch_size)
    error('sel2moon:BatchSize', 'batch_size must be a positive integer.');
end

% Scientific fingerprints exclude worker count and checkpoint frequency.
sources = {'sel2moon.build_database', 'sel2moon.prepare_seeds', ...
    'sel2moon.compute_seed', 'sel2rpo.compute_halo', ...
    'richardson_seed_halo_se3bp', 'richardson_coefficients_cr3bp', ...
    'libration_point_collinear_cr3bp', 'diffcorr_halo_halfperiod_se3bp', ...
    'rhs_SE3BP_3D', 'rhs_SE3BP_3D_STM', 'rhs_BCR4BP_3D', ...
    'se_cr3bp_to_emb_bcr4bp_state', 'emb_bcr4bp_to_se_cr3bp_state', 'rot2in'};
hashes = cellfun(@(name) rpo2ds.sha256_file(which(name)), sources, 'UniformOutput', false);
compatibility = struct('schema_version', 1, 'settings', settings, ...
    'sources', {sources}, 'source_hashes', {hashes});
completed = 0;
row_batches = {};
diagnostic_batches = {};
previous_seconds = 0;
timer = tic;
if isfile(checkpoint_file)
    saved = load(checkpoint_file, 'checkpoint');
    checkpoint = saved.checkpoint;
    if ~isequaln(checkpoint.compatibility, compatibility)
        error('sel2moon:CheckpointMismatch', ...
            ['Checkpoint settings or scientific code do not match:\n%s\n' ...
             'Explicitly discard this checkpoint to start from scratch.'], checkpoint_file);
    end
    context = checkpoint.context;
    completed = checkpoint.completed;
    row_batches = checkpoint.row_batches;
    diagnostic_batches = checkpoint.diagnostic_batches;
    previous_seconds = checkpoint.elapsed_seconds;
    fprintf('Resuming compatible checkpoint after %d seeds:\n%s\n', completed, checkpoint_file);
else
    context = sel2moon.prepare_seeds(settings);
end

% Save preparation too: a resumed campaign reuses exactly the same halo/STM.
if ~isfolder(fileparts(checkpoint_file)), mkdir(fileparts(checkpoint_file)); end
if ~isfolder(fileparts(database_file)), mkdir(fileparts(database_file)); end
checkpoint = struct('compatibility', compatibility, 'context', context, ...
    'completed', completed, 'row_batches', {row_batches}, ...
    'diagnostic_batches', {diagnostic_batches}, ...
    'elapsed_seconds', previous_seconds+toc(timer));
save_checkpoint(checkpoint_file, checkpoint);
use_parallel = rpo.parallel_available(execution.use_parallel);
if execution.use_parallel && ~use_parallel
    fprintf('Parallel Computing Toolbox is unavailable; using slower serial computation.\n');
end
if use_parallel
    pool = gcp;
    fprintf('SEL2-Moon: %d seeds; %d parallel workers.\n', settings.n_seed, pool.NumWorkers);
else
    fprintf('SEL2-Moon: %d seeds; serial computation.\n', settings.n_seed);
end

% Commit only complete batches; interruption reruns the unfinished batch.
while completed < settings.n_seed
    indices = completed+1:min(completed+execution.batch_size, settings.n_seed);
    rows = cell(numel(indices), 1);
    diagnostics = cell(numel(indices), 1);
    if use_parallel
        parfor index = 1:numel(indices)
            [rows{index}, diagnostics{index}] = sel2moon.compute_seed(indices(index), context);
        end
    else
        for index = 1:numel(indices)
            [rows{index}, diagnostics{index}] = sel2moon.compute_seed(indices(index), context);
        end
    end
    row_batches{end+1, 1} = vertcat(rows{:}); %#ok<AGROW>
    diagnostic_batches{end+1, 1} = struct2table(vertcat(diagnostics{:})); %#ok<AGROW>
    completed = indices(end);
    checkpoint.completed = completed;
    checkpoint.row_batches = row_batches;
    checkpoint.diagnostic_batches = diagnostic_batches;
    checkpoint.elapsed_seconds = previous_seconds+toc(timer);
    save_checkpoint(checkpoint_file, checkpoint);
    encounter_count = sum(cellfun(@numel, row_batches));
    fprintf('SEL2-Moon: %d/%d seeds (%.1f%%), %d encounters; %.1f s\n', ...
        completed, settings.n_seed, 100*completed/settings.n_seed, ...
        encounter_count, checkpoint.elapsed_seconds);
end

% Assign global indices after ordered concatenation, independent of batching.
table = struct2table(vertcat(row_batches{:}));
table.encounter_global_idx = (1:height(table)).';
per_seed = vertcat(diagnostic_batches{:});
diagnostics = struct('n_seeds', settings.n_seed, 'n_encounters', height(table), ...
    'n_hit_seeds', nnz(per_seed.n_encounters), ...
    'stage1_early_stop_count', nnz(per_seed.stage1_ended_early), ...
    'per_seed', per_seed, 'runtime_seconds', previous_seconds+toc(timer));
result_explorer = struct('settings', settings, 'orbit', context.orbit, ...
    'unstable', context.unstable, 'SEL2MoonDatabase', table, 'diagnostics', diagnostics);
cfg = struct('meta', struct('name', 'sel2_to_moon_database'));
result_explorer = rpo.finalize_result(result_explorer, cfg, ...
    'data_builders/build_sel2_to_moon_database.m', ...
    settings.orientation_input, 'generated_database');
result_explorer.provenance.scientific_fingerprint = compatibility;

temporary = [tempname(fileparts(database_file)), '.mat'];
save(temporary, 'result_explorer', '-v7.3');
movefile(temporary, database_file, 'f');
delete(checkpoint_file);
fprintf('Saved %d encounters from %d hit seeds to:\n%s\n', ...
    height(table), diagnostics.n_hit_seeds, database_file);
end

function save_checkpoint(filename, checkpoint)
%SAVE_CHECKPOINT Preserve completed batches until their replacement saves.
temporary = [tempname(fileparts(filename)), '.mat'];
save(temporary, 'checkpoint', '-v7');
movefile(temporary, filename, 'f');
end
