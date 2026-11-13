function selection = RPO2DS_time_capped_seed(me, rpo, idxME, I, J, nSol, time_policy)
%RPO2DS_TIME_CAPPED_SEED Rank paired seeds using the transfer-time lower bound.
%
% PROTOTYPE:
%   selection = RPO2DS_time_capped_seed(me, rpo, idxME, I, J, nSol, time_policy)
%
% DESCRIPTION:
%   Order the feasible pair pool monotonically by flight-time lower bound
%   and perigee cost. Prefer under-cap rows; retain the existing fallback
%   ranking when none is under cap. Legacy fields named tofExactDays retain
%   the screening value here; no exact flyby reconstruction is performed.
%
% INPUT:
%   me, rpo       Extracted Moon-Earth/RPO-Moon arrays from the patch matcher.
%   idxME         Maps local Moon-Earth neighbor indices to source rows.
%   I, J          Equal-length local pair indices into idxME and rpo.
%   nSol          Number of selected pairs.
%   time_policy   Struct containing max_tof_days and active time policy.
%
% OUTPUT:
%   selection     Struct with selected source/local indices, delta-v [km/s],
%                 time estimates [days], pool roles, ranking metadata, and
%                 reusable summary/refinement tables.

time_policy = local_normalize_time_policy(time_policy);

% Expand the feasible pair arrays into direct Moon-Earth and RPO row indices.
pair_me_global = idxME(double(I));
pair_rpo_local = double(J);
n_pairs = numel(pair_me_global);

% Read the pair-level Earth-flyby delta-v and the safe lower-bound transfer time.
pair_dvp_mag = me.dvpMag(pair_me_global);
pair_tof_lb_days = local_compute_tof_lower_bound_days(me.tof(pair_me_global), rpo.t_event(pair_rpo_local), time_policy);
pair_lb_keep = pair_tof_lb_days <= time_policy.max_tof_days;

% Build one fixed global pair ordering so the capped pool grows monotonically with Tmax.
master_pair_idx = local_rank_pairs_fixed(pair_tof_lb_days, pair_dvp_mag);

% Initialize one compact summary table reused by logs, tests, and driver scripts.
summary_table = table( ...
    (1:n_pairs).', ...
    pair_me_global(:), ...
    pair_rpo_local(:), ...
    rpo.row_index(pair_rpo_local), ...
    pair_dvp_mag(:), ...
    pair_tof_lb_days(:), ...
    pair_lb_keep(:), ...
    pair_tof_lb_days(:), ...
    NaN(n_pairs, 1), ...
    pair_lb_keep(:), ...
    false(n_pairs, 1), ...
    repmat(string("none"), n_pairs, 1), ...
    'VariableNames', {'candidateIndex', 'moonEarthIndex', 'rpoIndexLocal', 'rpoRowIndex', ...
    'deltavpMag', 'tofLowerBoundDays', 'tofLowerBoundKeep', 'tofExactDays', ...
    'flybyInboundDays', 'tofKeep', 'selected', 'poolRole'});

% Keep all lower-bound under-cap rows in the primary pool because no exact rebuild is performed anymore.
primary_pair_idx = master_pair_idx(pair_tof_lb_days(master_pair_idx) <= time_policy.max_tof_days);
summary_table.poolRole = local_mark_pool_role(summary_table.poolRole, primary_pair_idx, "primary");

% Select the best lower-bound under-cap candidates when they exist inside the full primary pool.
selected_pair_idx = zeros(0, 1);
if ~isempty(primary_pair_idx)
    selected_pair_idx = local_pick_selected_pairs(primary_pair_idx, summary_table, nSol, "under_cap");
end

% Fall back to the full feasible candidate set only when no lower-bound under-cap candidate exists.
fallback_pair_idx = zeros(0, 1);
selection_pool_role = "primary";
if isempty(selected_pair_idx)
    fallback_pair_idx = master_pair_idx;
    summary_table.poolRole = local_mark_pool_role(summary_table.poolRole, fallback_pair_idx, "fallback");
    selected_pair_idx = local_pick_selected_pairs(fallback_pair_idx, summary_table, nSol, "over_cap");
    selection_pool_role = "fallback";
end

% Flag the selected rows in the reusable screening summary.
summary_table.selected(selected_pair_idx) = true;

% Keep the full candidate pool visible in the same monotone lower-bound ordering.
refinement_pool_table = summary_table(master_pair_idx, :);

% Summarize the lower-bound screening pools and the final selection source compactly.
selected_pool_role = "none";
selected_candidate_count = 0;
if ~isempty(selected_pair_idx)
    selected_pool_role = selection_pool_role;
    if selection_pool_role == "primary"
        selected_candidate_count = numel(primary_pair_idx);
    elseif selection_pool_role == "fallback"
        selected_candidate_count = numel(fallback_pair_idx);
    end
end
pool_summary = table( ...
    n_pairs, ...
    numel(unique(pair_me_global)), ...
    numel(unique(pair_rpo_local)), ...
    numel(primary_pair_idx), ...
    numel(fallback_pair_idx), ...
    numel(primary_pair_idx), ...
    "tof_lower_bound_then_dvp", ...
    selected_pool_role, ...
    selected_candidate_count, ...
    numel(selected_pair_idx), ...
    'VariableNames', {'totalFeasiblePairs','totalFeasibleMoonEarthRows','totalFeasibleRPORows', ...
    'primaryPoolCount','fallbackPoolCount','primaryExactKeepCount', ...
    'poolOrdering','selectedPoolRole','selectedFromCandidateCount','selectedCount'});

% Pack the selected pair arrays used by the caller to build the patched seed table.
selection = struct();
selection.summary_table = summary_table;
selection.refinement_pool_table = refinement_pool_table;
selection.pool_summary = pool_summary;
selection.selected_pair_idx = selected_pair_idx(:);
selection.selected_me_global = pair_me_global(selected_pair_idx);
selection.selected_rpo_local = pair_rpo_local(selected_pair_idx);
selection.selected_deltavp_mag = pair_dvp_mag(selected_pair_idx);
selection.selected_tof_lower_bound_days = summary_table.tofLowerBoundDays(selected_pair_idx);
selection.selected_tof_exact_days = summary_table.tofExactDays(selected_pair_idx);
selection.selected_tof_keep = summary_table.tofKeep(selected_pair_idx);
selection.selected_flyby_inbound_days = summary_table.flybyInboundDays(selected_pair_idx);
selection.time_mode = char(time_policy.mode);
selection.max_tof_days = time_policy.max_tof_days;
selection.primary_pair_idx = primary_pair_idx(:);
selection.fallback_pair_idx = fallback_pair_idx(:);
selection.any_lower_bound_keep = any(pair_lb_keep);
selection.any_exact_keep = any(summary_table.tofKeep);
end

function pool_role = local_mark_pool_role(pool_role, pair_idx, role_label)
%LOCAL_MARK_POOL_ROLE Append the pool-role label without losing rows that belong to both pools.

for k = 1:numel(pair_idx)
    idx = pair_idx(k);
    if pool_role(idx) == "none"
        pool_role(idx) = role_label;
    elseif pool_role(idx) ~= role_label
        pool_role(idx) = "primary+fallback";
    end
end
end

function time_policy = local_normalize_time_policy(time_policy)
%LOCAL_NORMALIZE_TIME_POLICY Fill the capped-TOF policy defaults used by the seed selector.

if nargin < 1 || isempty(time_policy)
    time_policy = struct();
end
if ~isfield(time_policy, 'mode') || isempty(time_policy.mode)
    time_policy.mode = 'max_time';
end
if ~isfield(time_policy, 'max_tof_days') || isempty(time_policy.max_tof_days)
    time_policy.max_tof_days = inf;
end

% Use the same Earth-Moon scaling already used by the transfer solver.
time_policy.mode = lower(string(time_policy.mode));
time_policy.mu_EM = 0.01215;
time_policy.DU = 384400;
time_policy.TU = 1 / 2.665e-6;
time_policy.VU = time_policy.DU / time_policy.TU;
time_policy.TU_days = time_policy.TU / 86400;
time_policy.rSOI = (time_policy.mu_EM / (1 - time_policy.mu_EM))^(2/5);
end

function pair_tof_lb_days = local_compute_tof_lower_bound_days(T_lambert_sec, T_rpo_tu, time_policy)
%LOCAL_COMPUTE_TOF_LOWER_BOUND_DAYS Compute a safe lower bound on the absolute seed transfer time in days.

% Combine the signed RPO and Lambert leg times already available in the patch tables.
base_signed_days = T_rpo_tu(:) * time_policy.TU_days + T_lambert_sec(:) / 86400;
pair_tof_lb_days = abs(base_signed_days);

% The inbound Moon-flyby arc can shorten only the negative signed branch by at most one TU.
neg_idx = base_signed_days < 0;
pair_tof_lb_days(neg_idx) = max(abs(base_signed_days(neg_idx)) - time_policy.TU_days, 0);
end

function ranked_idx = local_rank_pairs_fixed(pair_tof_lb_days, pair_dvp_mag)
%LOCAL_RANK_PAIRS_FIXED Rank all feasible pairs once using the monotone capped-TOF ordering.

pair_idx = (1:numel(pair_tof_lb_days)).';
rank_tbl = table(pair_idx, pair_tof_lb_days(:), pair_dvp_mag(:), ...
    'VariableNames', {'pairIdx', 'tofLowerBoundDays', 'deltavpMag'});
rank_tbl = sortrows(rank_tbl, {'tofLowerBoundDays', 'deltavpMag', 'pairIdx'}, {'ascend', 'ascend', 'ascend'});
ranked_idx = rank_tbl.pairIdx;
end

function selected_pair_idx = local_pick_selected_pairs(pair_idx, summary_table, nSol, selection_mode)
%LOCAL_PICK_SELECTED_PAIRS Select the final capped-TOF rows from one exact-refined candidate pool.

pair_idx = pair_idx(:);
if isempty(pair_idx)
    selected_pair_idx = zeros(0, 1);
    return
end

% Use delta-v first for under-cap rows and shortest lower-bound time otherwise.
if selection_mode == "under_cap"
    ranked = table(pair_idx, summary_table.deltavpMag(pair_idx), summary_table.tofExactDays(pair_idx), ...
        'VariableNames', {'pairIdx', 'deltavpMag', 'tofExactDays'});
    ranked = sortrows(ranked, {'deltavpMag', 'tofExactDays', 'pairIdx'}, {'ascend', 'ascend', 'ascend'});
else
    ranked = table(pair_idx, summary_table.tofExactDays(pair_idx), summary_table.deltavpMag(pair_idx), ...
        'VariableNames', {'pairIdx', 'tofExactDays', 'deltavpMag'});
    ranked = sortrows(ranked, {'tofExactDays', 'deltavpMag', 'pairIdx'}, {'ascend', 'ascend', 'ascend'});
end
selected_pair_idx = ranked.pairIdx(1:min(nSol, height(ranked)));
end
