function phase_policy = RPO2DS_phase_policy(RPOMoonDatabase, params)
%RPO2DS_PHASE_POLICY Build the RPO database subset and wrapped phase metadata.
%
% PROTOTYPE:
%   phase_policy = RPO2DS_phase_policy(RPOMoonDatabase, params)
%
% DESCRIPTION:
%   Infer the periodic phase interval from orbitPhase samples. Preserve all
%   rows unless other_phase fixes the requested phase, in which case retain
%   the nearest two phase slices. Legacy tables without orbitPhase remain
%   on the fixed-phase branch and cannot satisfy other_phase requests.
%
% INPUT:
%   RPOMoonDatabase  RPO-Moon table; orbitPhase, when present, is in TU.
%   params           Struct with phase_mode (seed_phase/no_phase/other_phase)
%                      and optional phase_value [TU].
%
% OUTPUT:
%   phase_policy     Struct with database (including original row indices),
%                      phase grid/origin/step/period [TU], wrapped request,
%                      shortlist rows/distances, and phase-variable flags.

phase_mode = string(params.phase_mode);
if ~ismember('orbitPhase', RPOMoonDatabase.Properties.VariableNames)
    if phase_mode == "other_phase"
        error('RPO2DS:MissingPhaseData', ...
            'The selected database has no orbitPhase column for an other_phase request.');
    end
    % The legacy dense database has no phase samples. Keep its rows intact
    % and use the original fixed-phase solver branch.
    RPOMoonDatabasePrepared = RPOMoonDatabase;
    RPOMoonDatabasePrepared.rpoRowIndexOriginal = (1:height(RPOMoonDatabase)).';
    RPOMoonDatabasePrepared.phaseDistanceRequest = NaN(height(RPOMoonDatabase), 1);
    phase_policy = struct('mode', phase_mode, 'use_phase_variable', false, ...
        'phase_variable_mode', "explicit", 'database', RPOMoonDatabasePrepared, ...
        'phase_grid', 0, 'phase_origin', 0, 'phase_step', 2 * pi, ...
        'phase_period', 2 * pi, 'phase_request', [], ...
        'phase_request_wrapped', [], 'shortlist_rows', (1:height(RPOMoonDatabase)).', ...
        'shortlist_phases', zeros(height(RPOMoonDatabase), 1), ...
        'shortlist_distances', NaN(height(RPOMoonDatabase), 1), ...
        'phase_seed', 0, 'phase_seed_wrapped', 0);
    return
end
phase_values = RPOMoonDatabase.orbitPhase(:);
phase_grid = sort(unique(phase_values), 'ascend');

% Build the nominal phase step and the wrapped phase period from the database grid.
if numel(phase_grid) > 1
    phase_step = median(diff(phase_grid));
else
    phase_step = 2 * pi;
end
phase_origin = phase_grid(1);
phase_span = phase_grid(end) - phase_origin;
if local_has_periodic_endpoint(RPOMoonDatabase, phase_grid)
    phase_period = phase_span;
else
    phase_period = phase_span + phase_step;
end

% Wrap the requested phase onto the database interval when it is provided.
phase_request = [];
if isfield(params, 'phase_value') && ~isempty(params.phase_value)
    phase_request = local_wrap_phase(params.phase_value, phase_origin, phase_period);
end

% Precompute the wrapped distance from every database row to the requested phase.
phase_distance = NaN(size(phase_values));
if ~isempty(phase_request)
    phase_distance = local_wrapped_distance(phase_values, phase_request, phase_period);
end

% Start from the full database and attach the original row index.
RPOMoonDatabasePrepared = RPOMoonDatabase;
RPOMoonDatabasePrepared.rpoRowIndexOriginal = (1:height(RPOMoonDatabase)).'; % Keep the original row index for reports.
RPOMoonDatabasePrepared.phaseDistanceRequest = phase_distance; % Keep the wrapped distance to the requested phase.

% Restrict the database only when the user fixes a requested phase.
shortlist_rows = RPOMoonDatabasePrepared.rpoRowIndexOriginal;
if phase_mode == "other_phase"
    unique_phase_distance = local_wrapped_distance(phase_grid, phase_request, phase_period);
    [~, order_unique] = sortrows([unique_phase_distance(:), phase_grid(:)], [1 2]);
    shortlist_count = min(2, numel(order_unique)); % Keep the two closest phase grid values around the fixed NLP phase.
    shortlist_phase_values = phase_grid(order_unique(1:shortlist_count));
    keep_rows = false(height(RPOMoonDatabasePrepared), 1);
    for k = 1:numel(shortlist_phase_values)
        keep_rows = keep_rows | (abs(RPOMoonDatabasePrepared.orbitPhase - shortlist_phase_values(k)) <= 1e-12);
    end
    RPOMoonDatabasePrepared = RPOMoonDatabasePrepared(keep_rows, :);
    shortlist_rows = RPOMoonDatabasePrepared.rpoRowIndexOriginal;
end

% Pack the phase-policy output used by the transfer solver, tests, and script.
phase_policy = struct();
phase_policy.mode = phase_mode;
phase_policy.use_phase_variable = (phase_mode ~= "seed_phase");
phase_policy.phase_variable_mode = "explicit";
phase_policy.database = RPOMoonDatabasePrepared;
phase_policy.phase_grid = phase_grid;
phase_policy.phase_origin = phase_origin;
phase_policy.phase_step = phase_step;
phase_policy.phase_period = phase_period;
phase_policy.phase_request = phase_request;
phase_policy.phase_request_wrapped = phase_request;
phase_policy.shortlist_rows = shortlist_rows(:);
phase_policy.shortlist_phases = RPOMoonDatabasePrepared.orbitPhase(:);
phase_policy.shortlist_distances = RPOMoonDatabasePrepared.phaseDistanceRequest(:);
end

function phase_wrapped = local_wrap_phase(phase_raw, phase_origin, phase_period)
%LOCAL_WRAP_PHASE Wrap one phase value onto the database interval.

% Shift the phase by an integer number of wrapped periods.
phase_wrapped = phase_origin + mod(phase_raw - phase_origin, phase_period);
end

function dist = local_wrapped_distance(phase_values, phase_target, phase_period)
%LOCAL_WRAPPED_DISTANCE Compute the wrapped circular distance to one target phase.

% Keep the shortest signed distance around the wrapped phase interval.
delta = mod(phase_values - phase_target + 0.5 * phase_period, phase_period) - 0.5 * phase_period;
dist = abs(delta);
end

function has_endpoint = local_has_periodic_endpoint(RPOMoonDatabase, phase_grid)
%LOCAL_HAS_PERIODIC_ENDPOINT Detect databases whose final phase is the repeated periodic endpoint.

% A single phase value cannot expose whether the endpoint is repeated.
has_endpoint = false;
if numel(phase_grid) < 2
    return
end

% The RPO-Moon database is generated on these seed-design columns.
key_columns = {'deltaV', 'thetaV', 'verticalThetaV', 'n_hits'};
if ~all(ismember(key_columns, RPOMoonDatabase.Properties.VariableNames))
    return
end

% Compare the first and last phase slices after deterministic key sorting.
first_rows = RPOMoonDatabase(abs(RPOMoonDatabase.orbitPhase - phase_grid(1)) <= 1e-12, key_columns);
last_rows = RPOMoonDatabase(abs(RPOMoonDatabase.orbitPhase - phase_grid(end)) <= 1e-12, key_columns);
if height(first_rows) ~= height(last_rows)
    return
end
first_rows = sortrows(first_rows, key_columns);
last_rows = sortrows(last_rows, key_columns);

% Equal seed-design slices mean the last phase is the periodic endpoint.
has_endpoint = true;
for k = 1:numel(key_columns)
    column_name = key_columns{k};
    has_endpoint = has_endpoint && all(abs(first_rows.(column_name) - last_rows.(column_name)) <= 1e-12);
end
end
