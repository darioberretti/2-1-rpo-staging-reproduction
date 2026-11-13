function [RPOMoonDatabase, metadata] = build_database(settings, execution, database_file, checkpoint_file)
%BUILD_DATABASE Build an ordered encounter table with safe batch checkpoints.
%
% PROTOTYPE:
%   [RPOMoonDatabase, metadata] = rpomoon.build_database(settings, execution, ...
%       database_file, checkpoint_file)
%
% DESCRIPTION:
%   Enumerate the frozen impulse/phase grid, preserving repeated endpoints.
%   Resume only compatible scientific progress. Save a completed generated
%   database before deleting its checkpoint; never write reference data.
%
% INPUT:
%   settings         Resolved scientific settings/seed; schema and units are
%                    defined in the public scripts under data_builders/.
%   execution        use_parallel, batch_size, overwrite and producer_script
%                    (the producing script's repository-relative path).
%   database_file    Generated MAT destination, containing table plus metadata.
%   checkpoint_file  Builder-specific incomplete-progress MAT path.
%
% OUTPUT:
%   RPOMoonDatabase  [H x 20] Historical encounter schema, with H event rows.
%   metadata         Settings, units, source identities, counts and runtime.

% Protect the previous completed product before starting any lengthy work.
if isfile(database_file) && ~execution.overwrite
    if ~usejava('desktop')
        error('rpomoon:DatabaseExists', ...
            'Generated database exists: %s\nSet overwrite_generated_database = true to replace it.', ...
            database_file);
    end
    fprintf('The generated database will be replaced after successful completion:\n%s\n', database_file);
    answer = lower(strtrim(input('Proceed? [y/N]: ', 's')));
    if ~ismember(answer, {'y', 'yes'})
        error('rpomoon:Cancelled', 'Database replacement cancelled.');
    end
end

% Fingerprint the scientific implementation, independently of worker count.
sources = {'rpomoon.build_database', 'rpomoon.compute_task', 'propagate', ...
    'rhs_CR3BP_3D', 'rot2in'};
source_hashes = cellfun(@(name) rpo2ds.sha256_file(which(name)), sources, 'UniformOutput', false);
compatibility = struct('schema_version', 1, 'settings', settings, ...
    'sources', {sources}, 'source_hashes', {source_hashes});
row_batches = {};
completed = 0;
retry_count = 0;
previous_seconds = 0;
if isfile(checkpoint_file)
    saved = load(checkpoint_file, 'checkpoint');
    checkpoint = saved.checkpoint;
    if ~isequaln(checkpoint.compatibility, compatibility)
        error('rpomoon:CheckpointMismatch', ...
            ['Checkpoint settings, seed, schema or scientific code do not match:\n%s\n' ...
             'Explicitly discard this checkpoint, then rerun the builder to start from scratch.'], ...
            checkpoint_file);
    end
    row_batches = checkpoint.row_batches;
    completed = checkpoint.completed;
    retry_count = checkpoint.retry_count;
    previous_seconds = checkpoint.elapsed_seconds;
    fprintf('Resuming compatible checkpoint after %d completed tasks:\n%s\n', completed, checkpoint_file);
end

% Integer grid indices distinguish physically duplicated endpoint entries.
[magnitude, azimuth, elevation, phase] = ndgrid(settings.impulses_vu, ...
    settings.azimuths_rad, settings.elevations_rad, settings.phases_tu);
tasks = [magnitude(:), azimuth(:), elevation(:), phase(:)];
task_count = size(tasks, 1);
phase_seed = settings.seed_state_2d(:);
seed = [phase_seed(1:2); 0; phase_seed(3:4); 0];
departure_states = zeros(6, numel(settings.phases_tu));
for index = 1:numel(settings.phases_tu)
    [state, ~, ~, ~, ~, reached_phase] = propagate(seed, 0, settings.phases_tu(index), ...
        settings.mu, '3D', 'CR3BP', 'false');
    if abs(reached_phase-settings.phases_tu(index)) > 1e-12 || any(~isfinite(state))
        error('rpomoon:DeparturePropagation', 'RPO propagation did not reach phase %.17g TU.', ...
            settings.phases_tu(index));
    end
    departure_states(:, index) = state;
end
tasks_per_phase = numel(settings.impulses_vu)*numel(settings.azimuths_rad) ...
    * numel(settings.elevations_rad);
phase_indices = floor((0:task_count-1)/tasks_per_phase) + 1;
use_parallel = rpo.parallel_available(execution.use_parallel);
if execution.use_parallel && ~use_parallel
    fprintf('Parallel Computing Toolbox is unavailable; using slower serial computation.\n');
end
if use_parallel
    pool = gcp;
    fprintf('Encounter grid: %d tasks; %d parallel workers.\n', task_count, pool.NumWorkers);
else
    fprintf('Encounter grid: %d tasks; serial computation.\n', task_count);
end
mkdir_if_needed(fileparts(checkpoint_file));
mkdir_if_needed(fileparts(database_file));
timer = tic;
if ~isfile(checkpoint_file)
    checkpoint = struct('compatibility', compatibility, 'completed', completed, ...
        'row_batches', {row_batches}, 'retry_count', retry_count, 'elapsed_seconds', previous_seconds);
    save_checkpoint(checkpoint_file, checkpoint);
end

% Commit only completed batches. Failed/interrupted batches are rerun on resume.
while completed < task_count
    indices = completed+1:min(completed+execution.batch_size, task_count);
    rows = cell(numel(indices), 1);
    retried = false(numel(indices), 1);
    batch_tasks = tasks(indices, :);
    batch_states = departure_states(:, phase_indices(indices));
    if use_parallel
        parfor index = 1:numel(indices)
            [rows{index}, retried(index)] = rpomoon.compute_task( ...
                batch_tasks(index, :), batch_states(:, index), settings);
        end
    else
        for index = 1:numel(indices)
            [rows{index}, retried(index)] = rpomoon.compute_task( ...
                batch_tasks(index, :), batch_states(:, index), settings);
        end
    end
    row_batches{end+1, 1} = vertcat(rows{:}); %#ok<AGROW>
    completed = indices(end);
    retry_count = retry_count + nnz(retried);
    checkpoint = struct('compatibility', compatibility, 'completed', completed, ...
        'row_batches', {row_batches}, 'retry_count', retry_count, ...
        'elapsed_seconds', previous_seconds+toc(timer));
    save_checkpoint(checkpoint_file, checkpoint);
    encounter_count = sum(cellfun(@height, row_batches));
    fprintf('RPO-Moon: %d/%d tasks (%.1f%%), %d encounters, %d retries; %.1f s\n', ...
        completed, task_count, 100*completed/task_count, encounter_count, ...
        retry_count, previous_seconds+toc(timer));
end

% Save the compact completed product, then remove only this builder's cache.
RPOMoonDatabase = vertcat(row_batches{:});
metadata = struct('settings', settings, 'task_count', task_count, ...
    'encounter_count', height(RPOMoonDatabase), 'retry_count', retry_count, ...
    'runtime_seconds', previous_seconds+toc(timer), 'source_hashes', compatibility, ...
    'units', struct('impulse', 'VU', 'rotating_velocity', 'VU', ...
        'local_velocity', 'km/s', 'event_time', 'TU', 'event_time_days', 'days', ...
        'phase', 'TU', 'angles', 'rad except thetaMoon_deg [deg]'));
cfg = struct('meta', struct('name', 'rpo_moon_database'));
metadata = rpo.finalize_result(metadata, cfg, ...
    execution.producer_script, settings.seed_file, 'generated_database');
metadata.provenance.seed_sha256 = settings.seed_sha256;
temporary = [tempname(fileparts(database_file)), '.mat'];
save(temporary, 'RPOMoonDatabase', 'metadata', '-v7.3');
movefile(temporary, database_file, 'f');
if isfile(checkpoint_file)
    delete(checkpoint_file);
end
fprintf('Saved %d encounters to:\n%s\n', height(RPOMoonDatabase), database_file);
end

function mkdir_if_needed(directory)
%MKDIR_IF_NEEDED Create only the selected generated/cache parent directory.
if ~isfolder(directory)
    mkdir(directory);
end
end

function save_checkpoint(filename, checkpoint)
%SAVE_CHECKPOINT Preserve the previous valid batch until its replacement saves.
temporary = [tempname(fileparts(filename)), '.mat'];
save(temporary, 'checkpoint', '-v7.3');
movefile(temporary, filename, 'f');
end
