function [FullSolution, hMin, hBest, diagnostics, phase_policy, operational] = RPO2DS_operational_guess(vInfPlusEarth, RPOMoonDatabase, params, phase_template, constants)
%RPO2DS_OPERATIONAL_GUESS Build a patched seed under the active operational timing law.
%
% PROTOTYPE:
%   [FullSolution, hMin, hBest, diagnostics, phase_policy, operational] = ...
%       RPO2DS_operational_guess(vInfPlusEarth, RPOMoonDatabase, params, ...
%       phase_template, constants)
%
% DESCRIPTION:
%   Build an epoch-aware patch pool for initial-state or fixed-arrival
%   operation, reconstruct exact transfer timing, and select a seed compatible
%   with the specified onboard phase/epoch. The standard Run02 validation
%   path does not activate these optional operational policies.
%
% INPUT:
%   vInfPlusEarth   Outgoing Earth-relative seed velocity [km/s].
%   RPOMoonDatabase  Phase-aware RPO-Moon table.
%   params           Normalized RPO2DS_transfer parameters; physical time
%                      caps use days and geometry/phase epochs use TU.
%   phase_template   Phase grid/origin/period from RPO2DS_phase_policy.
%   constants        Transfer scaling and primary constants from RPO2DS_transfer.
%
% OUTPUT:
%   FullSolution     Selected patch table with operational timing annotations.
%   hMin, hBest      Feasible/selected Earth-perigee altitudes [km].
%   diagnostics      Patch-screening and timing diagnostics.
%   phase_policy     Updated phase metadata for the selected seed.
%   operational      Compact operational pool and selected-seed metadata.

guess_params = struct();
guess_params.N = params.N;
guess_params.hStart = params.hStart;
guess_params.sweep_rad = params.sweep_rad;
guess_params.nSol = params.operational_seed_pool_size;
guess_params.feasibleLimit = params.feasibleLimit;
guess_params.nMoonOrbitSamples = params.nMoonOrbitSamples;
guess_params.time_mode = 'max_time';
guess_params.maxTOF = local_build_pool_transfer_cap_days(params, phase_template);
guess_params.TU_days_operational = params.TU_days_operational;
guess_params.target_frame = params.target_frame;
guess_params.vinf_target_rth_km_s = params.vinf_target_rth_km_s;
guess_params.theta_sun0_mode = params.theta_sun0_mode;
guess_params.theta_sun0_value = params.theta_sun0_value;
guess_params.theta_sun0_seed = params.theta_sun0_seed;
guess_params.ws_sun_em_nd = params.ws_sun_em_nd;
guess_params.initial_epoch_time = params.initial_epoch_time;
guess_params.initial_rpo_phase = params.initial_rpo_phase;
guess_params.maxTransferTOF = params.maxTransferTOF;
if isfield(params, 'progress_logs')
    guess_params.progress_logs = params.progress_logs;
end
if isfield(params, 'rth_seed_epoch_nodes')
    guess_params.rth_seed_epoch_nodes = params.rth_seed_epoch_nodes;
end
if isfield(params, 'rth_seed_target_iterations')
    guess_params.rth_seed_target_iterations = params.rth_seed_target_iterations;
end
if isfield(params, 'rth_seed_nominal_lambert_tof_days')
    guess_params.rth_seed_nominal_lambert_tof_days = params.rth_seed_nominal_lambert_tof_days;
end
if isfield(params, 'rth_seed_max_lambert_tof_days')
    guess_params.rth_seed_max_lambert_tof_days = params.rth_seed_max_lambert_tof_days;
end
if isfield(params, 'rth_seed_max_pairs_per_epoch_node')
    guess_params.rth_seed_max_pairs_per_epoch_node = params.rth_seed_max_pairs_per_epoch_node;
end

% Let fixed-arrival mode build the Moon-Earth database from the final epoch.
if strcmp(params.operational_mode, 'fixed_arrival')
    guess_params.nSol = inf;
    guess_params.fixed_arrival_epoch_time = params.final_epoch_time;
    guess_params.fixed_arrival_phase_tolerance = params.fixed_arrival_phase_tolerance;
    guess_params.initial_epoch_time = params.initial_epoch_time;
    guess_params.initial_rpo_phase = params.initial_rpo_phase;
    guess_params.phase_origin = phase_template.phase_origin;
    guess_params.phase_period = phase_template.phase_period;
    guess_params.phase_step = phase_template.phase_step;
end

% Build the patched pool with the active target frame.
if local_use_time_aware_rth_guess(params)
    [FullSolutionPool, hMin, hBest, diagnostics] = ...
        RPO_Moon_Earth_guess_rth_time_aware(RPOMoonDatabase, guess_params, phase_template, constants);
else
    [FullSolutionPool, hMin, hBest, diagnostics] = RPO_Moon_Earth_guess(vInfPlusEarth, RPOMoonDatabase, guess_params);
end

% Fixed-arrival mode already selected the best feasible pair inside the patcher.
if strcmp(params.operational_mode, 'fixed_arrival')
    FullSolution = FullSolutionPool;
    timing_data = local_compute_exact_transfer_data(FullSolution, constants, phase_template.phase_period);
    FullSolution.seedToFSignedDays = timing_data.transfer_signed_days;
    FullSolution.seedToFExactDays = timing_data.transfer_forward_days;
    FullSolution.seedFlybyInboundDays = timing_data.flyby_inbound_days;
    FullSolution.seedPeriodLift = timing_data.period_lift_count;
    FullSolution.seedToFCapDays = repmat(guess_params.maxTOF, height(FullSolution), 1);
    FullSolution.seedToFKeep = FullSolution.seedToFExactDays <= guess_params.maxTOF + 1e-12;
    FullSolution.T_RPOMoonOriginal = FullSolution.T_RPOMoon;
    FullSolution.T_RPOMoon = FullSolution.T_RPOMoon + timing_data.period_lift_count * phase_template.phase_period;
    [FullSolution, selected_idx] = local_select_fixed_arrival_exact_seed(FullSolution, guess_params.maxTOF, params, phase_template);
    phase_seed_unwrapped = params.initial_rpo_phase + (params.final_epoch_time - params.initial_epoch_time) - ...
        FullSolution.seedToFExactDays(1) / params.TU_days_operational;
    screening = diagnostics.best_patch_summary;
    phase_policy = local_build_operational_phase_policy(params, phase_template, FullSolution, phase_seed_unwrapped, "from_initial_time");
    operational = local_build_operational_metadata(selected_idx, phase_seed_unwrapped, screening);
    return
end

% Lift every signed patched transfer onto one forward-time CR3BP branch.
timing_data = local_compute_exact_transfer_data(FullSolutionPool, constants, phase_template.phase_period);
FullSolutionPool.seedToFSignedDays = timing_data.transfer_signed_days;
FullSolutionPool.seedToFExactDays = timing_data.transfer_forward_days;
FullSolutionPool.seedFlybyInboundDays = timing_data.flyby_inbound_days;
FullSolutionPool.seedPeriodLift = timing_data.period_lift_count;
FullSolutionPool.seedToFCapDays = repmat(guess_params.maxTOF, height(FullSolutionPool), 1);
FullSolutionPool.seedToFKeep = FullSolutionPool.seedToFExactDays <= guess_params.maxTOF + 1e-12;
FullSolutionPool.T_RPOMoonOriginal = FullSolutionPool.T_RPOMoon;
FullSolutionPool.T_RPOMoon = FullSolutionPool.T_RPOMoon + timing_data.period_lift_count * phase_template.phase_period;

% Select one final candidate and the unwrapped departure phase used by the NLP.
if strcmp(params.operational_mode, 'initial_state')
    if local_use_time_aware_rth_guess(params)
        [selected_idx, screening, phase_seed_unwrapped] = local_select_initial_state_rth_seed(FullSolutionPool, params, phase_template);
    else
        [selected_idx, screening, phase_seed_unwrapped] = local_select_initial_state_seed(FullSolutionPool, params, phase_template);
    end
else
    [selected_idx, screening, phase_seed_unwrapped] = local_select_fixed_arrival_seed(FullSolutionPool, params, phase_template);
end

% Keep only the selected patched row for the shared downstream pipeline.
FullSolution = FullSolutionPool(selected_idx, :);
diagnostics.operational_screening = screening;

% Build the phase-policy payload reused by the shared NLP builder.
phase_policy = local_build_operational_phase_policy(params, phase_template, FullSolutionPool, phase_seed_unwrapped, "from_initial_time");

% Pack the compact operational metadata exposed by the result struct.
operational = local_build_operational_metadata(selected_idx, phase_seed_unwrapped, screening);
end

function tf = local_use_time_aware_rth_guess(params)
%LOCAL_USE_TIME_AWARE_RTH_GUESS Use the time-aware RTH seed search when the final epoch is not fixed.

tf = strcmpi(string(params.target_frame), "earth_sun_rth") && ...
    strcmpi(string(params.operational_mode), "initial_state");
end

function phase_policy = local_build_operational_phase_policy(params, phase_template, FullSolutionPool, phase_seed_unwrapped, phase_variable_mode)
%LOCAL_BUILD_OPERATIONAL_PHASE_POLICY Build the operational phase-policy payload reused by the shared NLP builder.

phase_policy = phase_template;
phase_policy.mode = string(params.operational_mode);
phase_policy.use_phase_variable = true;
phase_policy.phase_request = [];
phase_policy.phase_request_wrapped = [];
phase_policy.phase_seed = phase_seed_unwrapped;
phase_policy.phase_seed_wrapped = local_wrap_phase(phase_seed_unwrapped, phase_template.phase_origin, phase_template.phase_period);
phase_policy.phase_seed_nlp = phase_seed_unwrapped;
phase_policy.phase_bounds_nlp = local_build_phase_bounds(params, phase_template);
phase_policy.phase_variable_mode = phase_variable_mode;
phase_policy.shortlist_rows = FullSolutionPool.rpoRowIndex(:);
phase_policy.shortlist_phases = FullSolutionPool.phase(:);
phase_policy.shortlist_distances = NaN(height(FullSolutionPool), 1);
end

function operational = local_build_operational_metadata(selected_idx, phase_seed_unwrapped, screening)
%LOCAL_BUILD_OPERATIONAL_METADATA Pack the compact operational metadata exposed by the result struct.

operational = struct();
operational.selected_row = selected_idx;
operational.phase_seed_unwrapped = phase_seed_unwrapped;
operational.screening = screening;
end

function [FullSolution, selected_idx] = local_select_fixed_arrival_exact_seed(FullSolutionPool, cap_days, params, phase_template)
%LOCAL_SELECT_FIXED_ARRIVAL_EXACT_SEED Select the lowest-delta-v fixed-arrival row after exact timing and phase checks.

transfer_nd = FullSolutionPool.seedToFExactDays(:) / params.TU_days_operational;
available_time_nd = params.final_epoch_time - params.initial_epoch_time;
departure_phase_target = params.initial_rpo_phase + available_time_nd - transfer_nd;
departure_phase_wrapped = local_wrap_phase(departure_phase_target, phase_template.phase_origin, phase_template.phase_period);
departure_phase_mismatch = local_wrapped_distance(FullSolutionPool.phase(:), departure_phase_wrapped, phase_template.phase_period);
phase_tolerance = local_fixed_arrival_phase_tolerance(params, phase_template);
phase_keep = departure_phase_mismatch <= phase_tolerance;
under_cap_keep = FullSolutionPool.seedToFExactDays(:) <= cap_days + 1e-12;
candidate_idx = find(phase_keep & under_cap_keep);
if isempty(candidate_idx)
    candidate_idx = find(phase_keep);
end
ranked = table(candidate_idx, FullSolutionPool.deltavpMag(candidate_idx), FullSolutionPool.seedToFExactDays(candidate_idx), ...
    departure_phase_mismatch(candidate_idx), departure_phase_target(candidate_idx), departure_phase_wrapped(candidate_idx), ...
    under_cap_keep(candidate_idx), ...
    'VariableNames', {'idx', 'deltavpMag', 'tofDays', 'phaseMismatch', 'departurePhaseTarget', 'departurePhaseWrapped', 'underCap'});
ranked = sortrows(ranked, {'underCap', 'deltavpMag', 'tofDays', 'phaseMismatch', 'idx'}, {'descend', 'ascend', 'ascend', 'ascend', 'ascend'});

selected_idx = ranked.idx(1);
FullSolution = FullSolutionPool(selected_idx, :);
end

function phase_tolerance = local_fixed_arrival_phase_tolerance(params, phase_template)
%LOCAL_FIXED_ARRIVAL_PHASE_TOLERANCE Read the fixed-arrival departure-phase tolerance used by the exact selector.

phase_tolerance = phase_template.phase_step + 1e-12;
if isfield(params, 'fixed_arrival_phase_tolerance') && ~isempty(params.fixed_arrival_phase_tolerance)
    phase_tolerance = params.fixed_arrival_phase_tolerance;
end
end

function cap_days = local_build_pool_transfer_cap_days(params, phase_template)
%LOCAL_BUILD_POOL_TRANSFER_CAP_DAYS Build the transfer-time cap used while generating the patched pool.

% Keep the operational patched pool compatible with the later timing law.
if strcmp(params.operational_mode, 'fixed_arrival')
    available_time_nd = params.final_epoch_time - params.initial_epoch_time;
    cap_days = min(params.maxTransferTOF, max(available_time_nd, 0) * params.TU_days_operational);
else
    cap_days = params.maxTransferTOF;
end

% Never let the pool builder screen with a negative or zero duration.
cap_days = max(cap_days, phase_template.phase_step * params.TU_days_operational);
end

function [selected_idx, screening, phase_seed_unwrapped] = local_select_initial_state_seed(FullSolutionPool, params, phase_template)
%LOCAL_SELECT_INITIAL_STATE_SEED Rank the patched pool for the known-initial-state mode.

n_candidates = height(FullSolutionPool);
phase_period = phase_template.phase_period;
phase_wrapped = FullSolutionPool.phase(:);
phase_seed_rows = NaN(0, 1);
transfer_days_rows = NaN(0, 1);
wait_period_rows = NaN(0, 1);
wait_time_rows = NaN(0, 1);
total_time_rows = NaN(0, 1);
phase_mismatch_rows = NaN(0, 1);
dvp_rows = NaN(0, 1);
pool_row_idx = NaN(0, 1);
rpo_row_idx = NaN(0, 1);
keep_rows = false(0, 1);
phase_tolerance = local_phase_match_tolerance(params);
wait_cap_nd = local_initial_state_wait_cap_nd(params);

% Expand each patched candidate over the admissible full waiting periods.
for k = 1:n_candidates
    transfer_days = FullSolutionPool.seedToFExactDays(k);
    base_phase = phase_wrapped(k);
    n_min = max(ceil((params.initial_rpo_phase - base_phase) / phase_period), 0);
    wait0_nd = base_phase + n_min * phase_period - params.initial_rpo_phase;
    wait0_days = wait0_nd * params.TU_days_operational;
    wait_cap_days = wait_cap_nd * params.TU_days_operational;
    total_cap_periods = floor((params.maxTOF - transfer_days - wait0_days) / (phase_period * params.TU_days_operational)) + n_min;
    wait_cap_periods = floor((wait_cap_days - wait0_days) / (phase_period * params.TU_days_operational)) + n_min;
    n_max = min(total_cap_periods, wait_cap_periods);
    n_max = max(n_max, n_min);

    % Keep every waiting-period replica until the total cap is exceeded.
    for n_wait = n_min:n_max
        phase_unwrapped = base_phase + n_wait * phase_period;
        wait_time_nd = phase_unwrapped - params.initial_rpo_phase;
        wait_time_days = wait_time_nd * params.TU_days_operational;
        total_time_days = transfer_days + wait_time_days;
        departure_epoch_nd = params.initial_epoch_time + wait_time_nd;
        flyby_epoch_nd = departure_epoch_nd + FullSolutionPool.T_RPOMoon(k) + ...
            FullSolutionPool.seedFlybyInboundDays(k) / params.TU_days_operational;
        phase_mismatch = local_wrapped_distance(flyby_epoch_nd, FullSolutionPool.trueAnomalyMoonFlyby(k), 2 * pi);
        phase_seed_rows(end + 1, 1) = phase_unwrapped; %#ok<AGROW>
        transfer_days_rows(end + 1, 1) = transfer_days; %#ok<AGROW>
        wait_period_rows(end + 1, 1) = n_wait; %#ok<AGROW>
        wait_time_rows(end + 1, 1) = wait_time_days; %#ok<AGROW>
        total_time_rows(end + 1, 1) = total_time_days; %#ok<AGROW>
        phase_mismatch_rows(end + 1, 1) = phase_mismatch; %#ok<AGROW>
        dvp_rows(end + 1, 1) = FullSolutionPool.deltavpMag(k); %#ok<AGROW>
        pool_row_idx(end + 1, 1) = k; %#ok<AGROW>
        rpo_row_idx(end + 1, 1) = FullSolutionPool.rpoRowIndex(k); %#ok<AGROW>
        keep_rows(end + 1, 1) = (transfer_days <= params.maxTransferTOF + 1e-12) && ...
            (total_time_days <= params.maxTOF + 1e-12) && ...
            (wait_time_nd <= wait_cap_nd + 1e-12) && ...
            (phase_mismatch <= phase_tolerance + 1e-12); %#ok<AGROW>
    end
end

% Build the reusable screening table kept in the returned diagnostics.
screening = table((1:numel(pool_row_idx)).', pool_row_idx, rpo_row_idx, phase_seed_rows, ...
    wait_period_rows, transfer_days_rows, wait_time_rows, total_time_rows, phase_mismatch_rows, dvp_rows, ...
    keep_rows, false(numel(pool_row_idx), 1), ...
    'VariableNames', {'candidateReplicaIndex','candidatePoolIndex','rpoRowIndex','departurePhaseUnwrapped', ...
    'fullPeriodsWaited','transferTofDays','waitTimeDays','totalTimeDays','moonPhaseMismatch','deltavpMag','keep','selected'});

% Rank phase-consistent rows first because wrong inertial epoch breaks the NLP seed.
if any(keep_rows)
    ranked = sortrows(screening(keep_rows, :), {'moonPhaseMismatch', 'deltavpMag', 'totalTimeDays', 'transferTofDays', 'candidateReplicaIndex'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend', 'ascend'});
else
    ranked = sortrows(screening, {'moonPhaseMismatch', 'totalTimeDays', 'deltavpMag', 'candidateReplicaIndex'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend'});
end

selected_replica = ranked.candidateReplicaIndex(1);
screening.selected(selected_replica) = true;
selected_idx = screening.candidatePoolIndex(selected_replica);
phase_seed_unwrapped = screening.departurePhaseUnwrapped(selected_replica);
end

function wait_cap_nd = local_initial_state_wait_cap_nd(params)
%LOCAL_INITIAL_STATE_WAIT_CAP_ND Use the same two-RPO-period wait cap enforced by the NLP.

wait_cap_nd = inf;
if isfield(params, 'ws_sun_em_nd') && ~isempty(params.ws_sun_em_nd)
    wait_cap_nd = 2 * (-pi / params.ws_sun_em_nd);
end
end

function [selected_idx, screening, phase_seed_unwrapped] = local_select_initial_state_rth_seed(FullSolutionPool, params, phase_template)
%LOCAL_SELECT_INITIAL_STATE_RTH_SEED Rank the time-aware RTH patched pool by patch quality after wait-cap checks.

n_candidates = height(FullSolutionPool);
phase_period = phase_template.phase_period;
phase_wrapped = FullSolutionPool.phase(:);
phase_seed_rows = NaN(0, 1);
transfer_days_rows = NaN(0, 1);
wait_period_rows = NaN(0, 1);
wait_time_rows = NaN(0, 1);
total_time_rows = NaN(0, 1);
phase_mismatch_rows = NaN(0, 1);
dvp_rows = NaN(0, 1);
pool_row_idx = NaN(0, 1);
rpo_row_idx = NaN(0, 1);
keep_rows = false(0, 1);
wait_cap_nd = local_initial_state_wait_cap_nd(params);

% Expand each candidate over the admissible onboard wait replicas.
for k = 1:n_candidates
    transfer_days = FullSolutionPool.seedToFExactDays(k);
    base_phase = phase_wrapped(k);
    n_min = max(ceil((params.initial_rpo_phase - base_phase) / phase_period), 0);
    wait0_nd = base_phase + n_min * phase_period - params.initial_rpo_phase;
    wait0_days = wait0_nd * params.TU_days_operational;
    wait_cap_days = wait_cap_nd * params.TU_days_operational;
    total_cap_periods = floor((params.maxTOF - transfer_days - wait0_days) / (phase_period * params.TU_days_operational)) + n_min;
    wait_cap_periods = floor((wait_cap_days - wait0_days) / (phase_period * params.TU_days_operational)) + n_min;
    n_max = min(total_cap_periods, wait_cap_periods);
    n_max = max(n_max, n_min);

    for n_wait = n_min:n_max
        phase_unwrapped = base_phase + n_wait * phase_period;
        wait_time_nd = phase_unwrapped - params.initial_rpo_phase;
        wait_time_days = wait_time_nd * params.TU_days_operational;
        total_time_days = transfer_days + wait_time_days;
        departure_epoch_nd = params.initial_epoch_time + wait_time_nd;
        flyby_epoch_nd = departure_epoch_nd + FullSolutionPool.T_RPOMoon(k) + ...
            FullSolutionPool.seedFlybyInboundDays(k) / params.TU_days_operational;
        phase_mismatch = local_wrapped_distance(flyby_epoch_nd, FullSolutionPool.trueAnomalyMoonFlyby(k), 2 * pi);

        phase_seed_rows(end + 1, 1) = phase_unwrapped; %#ok<AGROW>
        transfer_days_rows(end + 1, 1) = transfer_days; %#ok<AGROW>
        wait_period_rows(end + 1, 1) = n_wait; %#ok<AGROW>
        wait_time_rows(end + 1, 1) = wait_time_days; %#ok<AGROW>
        total_time_rows(end + 1, 1) = total_time_days; %#ok<AGROW>
        phase_mismatch_rows(end + 1, 1) = phase_mismatch; %#ok<AGROW>
        dvp_rows(end + 1, 1) = FullSolutionPool.deltavpMag(k); %#ok<AGROW>
        pool_row_idx(end + 1, 1) = k; %#ok<AGROW>
        rpo_row_idx(end + 1, 1) = FullSolutionPool.rpoRowIndex(k); %#ok<AGROW>
        keep_rows(end + 1, 1) = (transfer_days <= params.maxTransferTOF + 1e-12) && ...
            (total_time_days <= params.maxTOF + 1e-12) && ...
            (wait_time_nd <= wait_cap_nd + 1e-12); %#ok<AGROW>
    end
end

screening = table((1:numel(pool_row_idx)).', pool_row_idx, rpo_row_idx, phase_seed_rows, ...
    wait_period_rows, transfer_days_rows, wait_time_rows, total_time_rows, phase_mismatch_rows, dvp_rows, ...
    keep_rows, false(numel(pool_row_idx), 1), ...
    'VariableNames', {'candidateReplicaIndex','candidatePoolIndex','rpoRowIndex','departurePhaseUnwrapped', ...
    'fullPeriodsWaited','transferTofDays','waitTimeDays','totalTimeDays','moonPhaseMismatch','deltavpMag','keep','selected'});

% The time-aware patcher already assigned each row to the closest epoch node.
if any(keep_rows)
    ranked = sortrows(screening(keep_rows, :), {'deltavpMag', 'moonPhaseMismatch', 'totalTimeDays', 'candidateReplicaIndex'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend'});
else
    ranked = sortrows(screening, {'totalTimeDays', 'moonPhaseMismatch', 'deltavpMag', 'candidateReplicaIndex'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend'});
end

selected_replica = ranked.candidateReplicaIndex(1);
screening.selected(selected_replica) = true;
selected_idx = screening.candidatePoolIndex(selected_replica);
phase_seed_unwrapped = screening.departurePhaseUnwrapped(selected_replica);
end

function [selected_idx, screening, phase_seed_unwrapped] = local_select_fixed_arrival_seed(FullSolutionPool, params, phase_template)
%LOCAL_SELECT_FIXED_ARRIVAL_SEED Rank the patched pool for the fixed-arrival operational mode.

transfer_days = FullSolutionPool.seedToFExactDays(:);
transfer_nd = transfer_days / params.TU_days_operational;
available_time_nd = params.final_epoch_time - params.initial_epoch_time;
departure_phase_target_unwrapped = params.initial_rpo_phase + available_time_nd - transfer_nd;
departure_phase_target_wrapped = local_wrap_phase(departure_phase_target_unwrapped, phase_template.phase_origin, phase_template.phase_period);
phase_mismatch = local_wrapped_distance(FullSolutionPool.phase(:), departure_phase_target_wrapped, phase_template.phase_period);
departure_epoch_nd = params.final_epoch_time - transfer_nd;
flyby_epoch_nd = departure_epoch_nd + FullSolutionPool.T_RPOMoon(:) + ...
    FullSolutionPool.seedFlybyInboundDays(:) / params.TU_days_operational;
moon_phase_mismatch = local_wrapped_distance(flyby_epoch_nd, FullSolutionPool.trueAnomalyMoonFlyby(:), 2 * pi);
phase_tolerance = local_phase_match_tolerance(params);
keep = (transfer_days <= (available_time_nd * params.TU_days_operational + 1e-12)) & ...
    (moon_phase_mismatch <= phase_tolerance + 1e-12);

% Build the reusable screening table kept in the returned diagnostics.
screening = table((1:height(FullSolutionPool)).', FullSolutionPool.rpoRowIndex(:), ...
    departure_phase_target_unwrapped, departure_phase_target_wrapped, transfer_days, ...
    phase_mismatch, moon_phase_mismatch, FullSolutionPool.deltavpMag(:), keep, false(height(FullSolutionPool), 1), ...
    'VariableNames', {'candidatePoolIndex','rpoRowIndex','departurePhaseUnwrapped', ...
    'departurePhaseWrapped','transferTofDays','departurePhaseMismatch','moonPhaseMismatch','deltavpMag','keep','selected'});

% Prefer rows that fit both the available transfer window and the Moon inertial phase.
if any(keep)
    ranked = sortrows(screening(keep, :), {'moonPhaseMismatch', 'departurePhaseMismatch', 'deltavpMag', 'transferTofDays', 'candidatePoolIndex'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend', 'ascend'});
else
    ranked = sortrows(screening, {'moonPhaseMismatch', 'departurePhaseMismatch', 'transferTofDays', 'deltavpMag', 'candidatePoolIndex'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend', 'ascend'});
end

selected_idx = ranked.candidatePoolIndex(1);
screening.selected(selected_idx) = true;
phase_seed_unwrapped = screening.departurePhaseUnwrapped(selected_idx);
end

function bounds = local_build_phase_bounds(params, ~)
%LOCAL_BUILD_PHASE_BOUNDS Build the unwrapped departure-phase bounds used by the operational NLP.

% Keep the departure phase after the known initial onboard phase.
lb = params.initial_rpo_phase;
if strcmp(params.operational_mode, 'fixed_arrival')
    ub = params.initial_rpo_phase + max(params.final_epoch_time - params.initial_epoch_time, 0);
else
    ub = params.initial_rpo_phase + params.maxTOF / params.TU_days_operational;
end

% Keep the interval valid even when the available time collapses to zero.
ub = max(ub, lb);
bounds = [lb, ub];
end

function phase_wrapped = local_wrap_phase(phase_raw, phase_origin, phase_period)
%LOCAL_WRAP_PHASE Wrap one phase value onto the database interval.

phase_wrapped = phase_origin + mod(phase_raw - phase_origin, phase_period);
end

function dist = local_wrapped_distance(phase_values, phase_target, phase_period)
%LOCAL_WRAPPED_DISTANCE Compute the wrapped circular distance to one target phase.

delta = mod(phase_values - phase_target + 0.5 * phase_period, phase_period) - 0.5 * phase_period;
dist = abs(delta);
end

function timing_data = local_compute_exact_transfer_data(FullSolutionPool, constants, phase_period)
%LOCAL_COMPUTE_EXACT_TRANSFER_DATA Rebuild one forward-time transfer duration for every patched candidate.

n_candidates = height(FullSolutionPool);
transfer_signed_nd = NaN(n_candidates, 1);
transfer_forward_nd = NaN(n_candidates, 1);
flyby_inbound_nd = NaN(n_candidates, 1);
period_lift_count = zeros(n_candidates, 1);

% Rebuild the signed CR3BP transfer and lift it onto the first forward branch.
for k = 1:n_candidates
    [transfer_signed_nd(k), flyby_inbound_nd(k)] = local_compute_exact_transfer_nd(FullSolutionPool(k, :), constants);
    [transfer_forward_nd(k), period_lift_count(k)] = local_lift_transfer_time(transfer_signed_nd(k), phase_period);
end

% Pack the signed and lifted timing data reused by the selector and the seed.
timing_data = struct();
timing_data.transfer_signed_days = transfer_signed_nd * constants.TU_days;
timing_data.transfer_forward_days = transfer_forward_nd * constants.TU_days;
timing_data.flyby_inbound_days = flyby_inbound_nd * constants.TU_days;
timing_data.period_lift_count = period_lift_count;
end

function [transfer_nd, flyby_inbound_nd] = local_compute_exact_transfer_nd(FullSolutionRow, constants)
%LOCAL_COMPUTE_EXACT_TRANSFER_ND Rebuild one signed patched transfer duration without assembling the full seed history.

% Rebuild the Moon flyby periapsis state from the paired asymptotic velocities.
[Rp, Vp] = local_get_perilune_state(FullSolutionRow.vInfMinusMoon, FullSolutionRow.vInfPlusMoon);
RpECI = RTH2ECI_moon(FullSolutionRow.RI, Rp, [], 1);
VpECI = RTH2ECI_moon(FullSolutionRow.RI, Vp, [], 0);
x0 = [RpECI(:) / constants.DU; VpECI(:) / constants.VU];
time = FullSolutionRow.trueAnomalyMoonFlyby;
xx0_peri = in2rot(x0.', time, constants.mu_EM);

% Propagate only the inbound flyby branch because the exact transfer time
% needs the Moon-peri to SOI-entry duration used by the assembled seed.
opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
    'Events', @(t, x) local_soi_event_cr3bp(t, x, constants.rSOI, constants.mu_EM));
[tB, ~] = ode113(@(t, x) rhs_CR3BP_3D(t, x, constants.mu_EM), [time, time - 1], xx0_peri, opts_ode);
flyby_inbound_nd = tB(1) - tB(end);

% Keep the signed patched transfer before lifting it onto the forward branch.
transfer_nd = FullSolutionRow.T_RPOMoon + flyby_inbound_nd + FullSolutionRow.T_lambert / constants.TU;
end

function phase_tolerance = local_phase_match_tolerance(params)
%LOCAL_PHASE_MATCH_TOLERANCE Use the Moon sampling density as the admissible inertial phase mismatch.

phase_tolerance = pi / max(params.nMoonOrbitSamples, 1);
end

function [transfer_forward_nd, period_lift_count] = local_lift_transfer_time(transfer_signed_nd, phase_period)
%LOCAL_LIFT_TRANSFER_TIME Lift one signed transfer duration onto the smallest forward-time branch.

period_lift_count = 0;
transfer_forward_nd = transfer_signed_nd;

% Add whole RPO periods until the transfer becomes forward in time.
if transfer_forward_nd <= 0
    period_lift_count = ceil((-transfer_forward_nd + 1e-12) / phase_period);
    transfer_forward_nd = transfer_forward_nd + period_lift_count * phase_period;
end

% Keep a nonzero forward duration when the signed time lands on one wrapped boundary.
if transfer_forward_nd <= 1e-12
    period_lift_count = period_lift_count + 1;
    transfer_forward_nd = transfer_forward_nd + phase_period;
end
end

function [Rp, Vp] = local_get_perilune_state(v_in, v_out)
%LOCAL_GET_PERILUNE_STATE Rebuild the perilune state from the incoming and outgoing Moon-relative vInf vectors.

muM = astroConstants(20);
dotProd = dot(v_in, v_out);
norm2 = norm(v_in) * norm(v_out);
delta = acos(dotProd / norm2);
v_mid = (norm(v_in) + norm(v_out)) / 2;
RpMag = muM / v_mid^2 * (1 / sin(delta / 2) - 1);
VpMag = sqrt(v_mid^2 + 2 * muM / RpMag);
Rp = RpMag / (2 * v_mid * sin(delta / 2)) * (v_in - v_out);
Vp = VpMag / (2 * v_mid * cos(delta / 2)) * (v_in + v_out);
end

function [value, isterminal, direction] = local_soi_event_cr3bp(~, x, rSOI, mu)
%LOCAL_SOI_EVENT_CR3BP Stop the inbound flyby branch at the Moon SOI crossing.

posMoon = [1 - mu, 0, 0].';
value = norm(x(1:3) - posMoon) - rSOI;
isterminal = 1;
direction = 0;
end
