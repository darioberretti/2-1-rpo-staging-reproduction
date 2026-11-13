function [solutionTable, diagnostics] = MoonEarth_RPOMoon_cellVersion_opt_geometry(MoonEarthDatabase, RPOMoonDatabase, nSol, feasibleLimit, hMinMoonFlyby, timePolicy, geometryPolicy)
%MOONEARTH_RPOMOON_CELLVERSION_OPT_GEOMETRY Match Moon-Earth patches to the RPO-Moon database.
%
% PROTOTYPE:
%   [solutionTable, diagnostics] = MoonEarth_RPOMoon_cellVersion_opt_geometry( ...
%       MoonEarthDatabase, RPOMoonDatabase, nSol, feasibleLimit, ...
%       hMinMoonFlyby, timePolicy, geometryPolicy)
%
% DESCRIPTION:
%   Batch candidates by the prescribed initial geometry, then apply the
%   shared lunar speed and flyby-deflection filters. Rank surviving pairs
%   by the active time/perigee-delta-v policy and retain selected patch rows.
%
% INPUT:
%   MoonEarthDatabase  Table from interplanetary_MoonEarth_patch_opt or its
%                        epoch-aware variant; physical vectors use km/km/s.
%   RPOMoonDatabase   RPO-Moon table; phase/t_event use normalized TU and
%                        t_event_days uses days. See docs/DATA.md.
%   nSol              Requested number of selected rows.
%   feasibleLimit     Optional lunar speed-match tolerance [km/s], default 1e-3.
%   hMinMoonFlyby     Optional minimum lunar altitude [km], default 100.
%   timePolicy        Optional mode (default optimization_only) and
%                        max_tof_days (default Inf).
%   geometryPolicy    Optional prescribed_geometry mode and
%                        initial_geometry_time [TU], default 0.
%
% OUTPUT:
%   solutionTable     Selected patch candidates, or the standard empty table.
%   diagnostics       Prefilter/flyby counts and timing/selection summaries
%                        including per-phase batch counts.

if nargin < 5 || isempty(hMinMoonFlyby)
    hMinMoonFlyby = 100; % [km]
end
if nargin < 4 || isempty(feasibleLimit)
    feasibleLimit = 1e-3; % [km/s]
end
if nargin < 6 || isempty(timePolicy)
    timePolicy = struct();
end
if nargin < 7 || isempty(geometryPolicy)
    geometryPolicy = struct();
end
if ~isfield(timePolicy, 'mode') || isempty(timePolicy.mode)
    timePolicy.mode = 'optimization_only';
end
if ~isfield(timePolicy, 'max_tof_days') || isempty(timePolicy.max_tof_days)
    timePolicy.max_tof_days = inf;
end
if ~isfield(geometryPolicy, 'mode') || isempty(geometryPolicy.mode)
    geometryPolicy.mode = 'prescribed_geometry';
end
if ~isfield(geometryPolicy, 'initial_geometry_time') || isempty(geometryPolicy.initial_geometry_time)
    geometryPolicy.initial_geometry_time = 0;
end

% Load the Moon flyby constants once.
persistent MU_MOON R_MOON
if isempty(MU_MOON)
    MU_MOON = astroConstants(20);
    R_MOON  = astroConstants(30);
end
muM = MU_MOON;
rpmin = R_MOON + hMinMoonFlyby;

% Extract the Earth-side and RPO-side arrays used by the batched comparison.
me = local_extract_moon_earth(MoonEarthDatabase);
rpo = local_extract_rpo_moon(RPOMoonDatabase, geometryPolicy);
me_total = numel(me.vInfMag);
rpo_total = numel(rpo.vInfMag);

% Apply one global Moon-Earth |vInf| prefilter before any phase batching.
[me.keep_global, diagnostics] = local_prefilter_by_vinf(me.vInfMag, rpo.vInfMag, feasibleLimit);
if isempty(me.keep_global)
    batch_summary = local_build_phase_summary(0);
    solutionTable = local_empty_solution_table();
    diagnostics = local_build_empty_diagnostics(batch_summary, geometryPolicy, diagnostics.vInfMoon_within_constraints);
    diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(me.keep_global), 0, 0, 0, 0, 0, 0, 0);
    return
end

% Build the phase samples and prebuild the row batches once.
[phase_samples, me.phase_sample_idx] = local_build_phase_batches(me.theta);
rpo.phase_sample_idx = local_assign_phase_batches(rpo.theta_corrected, phase_samples);
[me.phase_batches, batch_summary] = local_build_phase_batch_cells(me.phase_sample_idx, numel(phase_samples), me.keep_global);
[rpo.phase_batches, ~] = local_build_phase_batch_cells(rpo.phase_sample_idx, numel(phase_samples), uint32((1:numel(rpo.vInfMag)).'));
batch_summary.phaseSampleValue = phase_samples(:);

if strcmpi(timePolicy.mode, 'max_time')
    % Keep the original capped-time seed selection when the capped-TOF mode is active.
    [pair_me_global, pair_rpo_local, batch_summary, filter_counts] = ...
        local_collect_batched_pairs(me, rpo, feasibleLimit, rpmin, muM, me.phase_batches, rpo.phase_batches, batch_summary);
    if isempty(pair_me_global)
        solutionTable = local_empty_solution_table();
        diagnostics = local_build_empty_diagnostics(batch_summary, geometryPolicy, diagnostics.vInfMoon_within_constraints);
        diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(me.keep_global), ...
            filter_counts.moonEarthRowsWithNeighbors, filter_counts.candidatePairsBeforeFlyby, filter_counts.candidatePairsAfterFlyby, ...
            filter_counts.uniqueMoonEarthRowsBeforeFlyby, filter_counts.uniqueRPORowsBeforeFlyby, ...
            filter_counts.uniqueMoonEarthRowsAfterFlyby, filter_counts.uniqueRPORowsAfterFlyby);
        return
    end
    diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(me.keep_global), ...
        filter_counts.moonEarthRowsWithNeighbors, filter_counts.candidatePairsBeforeFlyby, filter_counts.candidatePairsAfterFlyby, ...
        filter_counts.uniqueMoonEarthRowsBeforeFlyby, filter_counts.uniqueRPORowsBeforeFlyby, ...
        filter_counts.uniqueMoonEarthRowsAfterFlyby, filter_counts.uniqueRPORowsAfterFlyby);

    idxME = unique(pair_me_global, 'stable');
    [~, I_local] = ismember(pair_me_global, idxME);
    timeSelection = RPO2DS_time_capped_seed(me, rpo, idxME, uint32(I_local(:)), uint32(pair_rpo_local(:)), nSol, timePolicy);
    if isempty(timeSelection.selected_me_global)
        solutionTable = local_empty_solution_table();
        diagnostics = local_build_empty_diagnostics(batch_summary, geometryPolicy, diagnostics.vInfMoon_within_constraints);
        diagnostics.time_summary = timeSelection.summary_table;
        diagnostics.time_policy = rmfield(timeSelection, 'summary_table');
        diagnostics.refinement_pool_table = timeSelection.refinement_pool_table;
        diagnostics.refinement_pool_summary = timeSelection.pool_summary;
        diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(me.keep_global), ...
            filter_counts.moonEarthRowsWithNeighbors, filter_counts.candidatePairsBeforeFlyby, filter_counts.candidatePairsAfterFlyby, ...
            filter_counts.uniqueMoonEarthRowsBeforeFlyby, filter_counts.uniqueRPORowsBeforeFlyby, ...
            filter_counts.uniqueMoonEarthRowsAfterFlyby, filter_counts.uniqueRPORowsAfterFlyby);
        return
    end

    diagnostics.time_summary = timeSelection.summary_table;
    diagnostics.time_policy = rmfield(timeSelection, 'summary_table');
    diagnostics.refinement_pool_table = timeSelection.refinement_pool_table;
    diagnostics.refinement_pool_summary = timeSelection.pool_summary;
    solutionTable = local_pack_solution_table( ...
        me, rpo, timeSelection.selected_me_global, timeSelection.selected_rpo_local, ...
        timeSelection.selected_deltavp_mag, RPOMoonDatabase, diagnostics.time_policy, geometryPolicy, phase_samples);
else
    % Scan Moon-Earth rows in increasing Earth-flyby delta-v and stop as soon as the best nSol rows are feasible.
    [chosenME_global, RPOMoonBest_idx, minDeltavp, batch_summary, filter_counts] = ...
        local_select_best_pairs_by_dvp(me, rpo, feasibleLimit, rpmin, muM, nSol, me.phase_batches, rpo.phase_batches, batch_summary);
    if isempty(chosenME_global)
        solutionTable = local_empty_solution_table();
        diagnostics = local_build_empty_diagnostics(batch_summary, geometryPolicy, diagnostics.vInfMoon_within_constraints);
        diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(me.keep_global), ...
            filter_counts.moonEarthRowsWithNeighbors, filter_counts.candidatePairsBeforeFlyby, filter_counts.candidatePairsAfterFlyby, ...
            filter_counts.uniqueMoonEarthRowsBeforeFlyby, filter_counts.uniqueRPORowsBeforeFlyby, ...
            filter_counts.uniqueMoonEarthRowsAfterFlyby, filter_counts.uniqueRPORowsAfterFlyby);
        return
    end
    diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(me.keep_global), ...
        filter_counts.moonEarthRowsWithNeighbors, filter_counts.candidatePairsBeforeFlyby, filter_counts.candidatePairsAfterFlyby, ...
        filter_counts.uniqueMoonEarthRowsBeforeFlyby, filter_counts.uniqueRPORowsBeforeFlyby, ...
        filter_counts.uniqueMoonEarthRowsAfterFlyby, filter_counts.uniqueRPORowsAfterFlyby);

    % Reuse the legacy diagnostics schema for the non-capped mode.
    diagnostics.time_summary = table();
    diagnostics.time_policy = struct('time_mode', char(lower(string(timePolicy.mode))), ...
        'selected_tof_lower_bound_days', NaN(size(minDeltavp(:))), ...
        'selected_tof_exact_days', NaN(size(minDeltavp(:))), ...
        'selected_tof_keep', false(size(minDeltavp(:))), ...
        'selected_flyby_inbound_days', NaN(size(minDeltavp(:))), ...
        'max_tof_days', timePolicy.max_tof_days);
    diagnostics.refinement_pool_table = table();
    diagnostics.refinement_pool_summary = table();
    solutionTable = local_pack_solution_table(me, rpo, chosenME_global, RPOMoonBest_idx, minDeltavp, RPOMoonDatabase, diagnostics.time_policy, geometryPolicy, phase_samples);
end

% Store the geometry-aware screening diagnostics for tests and logs.
diagnostics.geometry_policy = geometryPolicy;
diagnostics.phase_samples = phase_samples(:);
diagnostics.phase_batch_summary = batch_summary;
end

function me = local_extract_moon_earth(MoonEarthDatabase)
%LOCAL_EXTRACT_MOON_EARTH Extract the Moon-Earth arrays once to avoid repeated table indexing.

me = struct();
me.vInfMag = MoonEarthDatabase.vInfMoonMag;
me.vInf    = MoonEarthDatabase.vInfMoon;
me.dvpvec  = MoonEarthDatabase.deltavp;
me.dvpMag  = MoonEarthDatabase.deltavpMag;
me.Vp      = MoonEarthDatabase.VF; % Keep the legacy packed variable used by the patch table.
me.R2      = MoonEarthDatabase.R2;
me.theta   = MoonEarthDatabase.thetaMoon;
me.tof     = MoonEarthDatabase.ToF;
me.RI      = MoonEarthDatabase.R1;
me.VI      = MoonEarthDatabase.VI;
me.moonEpochNd = local_optional_column(MoonEarthDatabase, "moonEpochNd");
me.earthEpochNd = local_optional_column(MoonEarthDatabase, "earthEpochNd");
me.epochNodeIndex = local_optional_column(MoonEarthDatabase, "epochNodeIndex");
me.earthPeriSampleIndex = local_optional_column(MoonEarthDatabase, "earthPeriSampleIndex");
end

function value = local_optional_column(tbl, columnName)
%LOCAL_OPTIONAL_COLUMN Read optional time-aware metadata without changing legacy ECI tables.

value = [];
if ismember(columnName, string(tbl.Properties.VariableNames))
    value = tbl.(columnName);
end
end

function rpo = local_extract_rpo_moon(RPOMoonDatabase, geometryPolicy)
%LOCAL_EXTRACT_RPO_MOON Extract the RPO-Moon arrays and add the corrected Moon phase used by the geometry-aware matcher.

rpo = struct();
rpo.vInfMag = RPOMoonDatabase.vInfMagMoon;
rpo.vInf    = RPOMoonDatabase.vInfMoon_RTH;
rpo.dv      = RPOMoonDatabase.deltaV;
rpo.dv0     = RPOMoonDatabase.deltaV_vect;
rpo.t_event = RPOMoonDatabase.t_event;
rpo.theta_raw = wrapTo2Pi(deg2rad(RPOMoonDatabase.thetaMoon_deg));
rpo.theta_corrected = wrapTo2Pi(geometryPolicy.initial_geometry_time + rpo.t_event);
rpo.row_index = (1:height(RPOMoonDatabase)).';
if ismember("rpoRowIndexOriginal", RPOMoonDatabase.Properties.VariableNames)
    rpo.row_index = RPOMoonDatabase.rpoRowIndexOriginal;
end
if ismember("phaseDistanceRequest", RPOMoonDatabase.Properties.VariableNames)
    rpo.phase_distance = RPOMoonDatabase.phaseDistanceRequest;
else
    rpo.phase_distance = NaN(height(RPOMoonDatabase), 1);
end
if ismember("orbitPhase", RPOMoonDatabase.Properties.VariableNames)
    rpo.phase = RPOMoonDatabase.orbitPhase;
else
    rpo.phase = zeros(height(RPOMoonDatabase), 1);
end
end

function [phase_samples, phase_idx] = local_build_phase_batches(theta_values)
%LOCAL_BUILD_PHASE_BATCHES Build the sorted Moon-Earth phase samples and assign each Moon-Earth row to its nearest sample.

theta_wrapped = wrapTo2Pi(theta_values(:));
phase_samples = unique(theta_wrapped, 'stable');
phase_samples = sort(phase_samples);
phase_idx = local_assign_phase_batches(theta_wrapped, phase_samples);
end

function phase_idx = local_assign_phase_batches(theta_values, phase_samples)
%LOCAL_ASSIGN_PHASE_BATCHES Assign each wrapped phase to the nearest Moon-Earth phase sample on the circle.

theta_wrapped = wrapTo2Pi(theta_values(:));
phase_samples = wrapTo2Pi(phase_samples(:).');
angle_diff = abs(angle(exp(1i * (theta_wrapped - phase_samples))));
[~, phase_idx] = min(angle_diff, [], 2);
phase_idx = uint32(phase_idx);
end

function [phase_batches, batch_summary] = local_build_phase_batch_cells(phase_idx, n_batches, row_idx)
%LOCAL_BUILD_PHASE_BATCH_CELLS Prebuild the row lists that belong to each phase sample.

phase_batches = cell(n_batches, 1);
batch_summary = local_build_phase_summary(n_batches);
row_idx = uint32(row_idx(:));
phase_idx = uint32(phase_idx(:));
for i = 1:n_batches
    phase_batches{i} = row_idx(phase_idx(double(row_idx)) == i);
end
end

function batch_summary = local_build_phase_summary(n_batches)
%LOCAL_BUILD_PHASE_SUMMARY Build the phase-batch diagnostics table.

batch_summary = table('Size', [n_batches, 11], ...
    'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'phaseSampleIndex', 'phaseSampleValue', 'moonEarthCount', 'rpoCount', ...
    'pairCountBeforeFlyby','pairCountAfterFlyby','pairCount', ...
    'moonEarthRowsBeforeFlyby','rpoRowsBeforeFlyby','moonEarthRowsAfterFlyby','rpoRowsAfterFlyby'});
batch_summary.phaseSampleIndex = (1:n_batches).';
batch_summary.phaseSampleValue = zeros(n_batches, 1);
end

function [pair_me_global, pair_rpo_local, batch_summary, filter_counts] = local_collect_batched_pairs(me, rpo, feasibleLimit, rpmin, muM, me_phase_batches, rpo_phase_batches, batch_summary)
%LOCAL_COLLECT_BATCHED_PAIRS Build the feasible pair list by comparing only rows that share the same phase batch.

n_batches = numel(me_phase_batches);
pair_me_global = zeros(0, 1, 'uint32');
pair_rpo_local = zeros(0, 1, 'uint32');
pair_me_before_global = zeros(0, 1, 'uint32');
pair_rpo_before_local = zeros(0, 1, 'uint32');

for i = 1:n_batches
    me_idx = me_phase_batches{i};
    rpo_idx = rpo_phase_batches{i};
    batch_summary.moonEarthCount(i) = numel(me_idx);
    batch_summary.rpoCount(i) = numel(rpo_idx);

    % Skip empty phase batches early.
    if isempty(me_idx) || isempty(rpo_idx)
        continue
    end

    % Search only inside the current phase batch for close |vInf| neighbors.
    me_local_mags = me.vInfMag(double(me_idx));
    neigh = local_find_neighbors_by_magnitude(rpo.vInfMag(rpo_idx), me_local_mags, feasibleLimit);
    counts = cellfun('length', neigh);
    has_any = counts > 0;
    if ~any(has_any)
        continue
    end

    idxME_local = uint32(find(has_any));
    neigh = neigh(has_any);
    counts = counts(has_any);

    % Flatten the batch-local pair arrays before the flyby-feasibility filter.
    I = repelem(uint32((1:numel(idxME_local))'), counts);
    J = vertcat(neigh{:});
    me_batch_global = uint32(me_idx(double(idxME_local(double(I)))));
    rpo_batch_local = uint32(rpo_idx(double(J)));
    pair_me_before_global = [pair_me_before_global; me_batch_global(:)]; %#ok<AGROW>
    pair_rpo_before_local = [pair_rpo_before_local; rpo_batch_local(:)]; %#ok<AGROW>

    % Apply the same flyby-deflection feasibility filter used by the legacy matcher.
    [vPlus, vMinus, magPlus, magMinus] = local_gather_candidate_vectors(me, rpo, me_batch_global, rpo_batch_local);
    ok = local_deflection_feasible(vPlus, vMinus, magPlus, magMinus, rpmin, muM);
    me_batch_global = me_batch_global(ok);
    rpo_batch_local = rpo_batch_local(ok);

    batch_summary.pairCountBeforeFlyby(i) = numel(I);
    batch_summary.pairCountAfterFlyby(i) = numel(me_batch_global);
    batch_summary.pairCount(i) = batch_summary.pairCountAfterFlyby(i);
    batch_summary.moonEarthRowsBeforeFlyby(i) = numel(unique(double(me_idx(double(idxME_local)))));
    batch_summary.rpoRowsBeforeFlyby(i) = numel(unique(double(rpo_idx(double(J)))));
    batch_summary.moonEarthRowsAfterFlyby(i) = numel(unique(double(me_batch_global)));
    batch_summary.rpoRowsAfterFlyby(i) = numel(unique(double(rpo_batch_local)));
    pair_me_global = [pair_me_global; me_batch_global(:)]; %#ok<AGROW>
    pair_rpo_local = [pair_rpo_local; rpo_batch_local(:)]; %#ok<AGROW>
end

filter_counts = local_build_filter_counts(pair_me_before_global, pair_rpo_before_local, pair_me_global, pair_rpo_local);
end

function [idxME_all, diagnostics] = local_prefilter_by_vinf(me_vInfMag, rpo_vInfMag, feasibleLimit)
%LOCAL_PREFILTER_BY_VINF Prefilter one Moon-Earth batch by the admissible RPO |vInf| band.

lo = min(rpo_vInfMag) - feasibleLimit;
hi = max(rpo_vInfMag) + feasibleLimit;
keepME = (me_vInfMag >= lo) & (me_vInfMag <= hi);
idxME_all = uint32(find(keepME));
diagnostics = struct('vInfMoon_within_constraints', numel(idxME_all));
end

function neigh = local_find_neighbors_by_magnitude(rpo_vInfMag, me_mags, feasibleLimit)
%LOCAL_FIND_NEIGHBORS_BY_MAGNITUDE Find |vInf| neighbors with the same fast KD-tree fallback logic used by the legacy matcher.

if exist('createns', 'file') == 2 && exist('rangesearch', 'file') == 2
    mdl = createns(rpo_vInfMag(:), 'NSMethod', 'kdtree');
    neigh = rangesearch(mdl, me_mags(:), feasibleLimit);
else
    [rpo_sorted, perm] = sort(rpo_vInfMag(:));
    neigh = cell(numel(me_mags), 1);
    j0 = 1;
    for ii = 1:numel(me_mags)
        x = me_mags(ii);
        while j0 <= numel(rpo_sorted) && rpo_sorted(j0) < x - feasibleLimit
            j0 = j0 + 1;
        end
        j = j0;
        while j <= numel(rpo_sorted) && rpo_sorted(j) <= x + feasibleLimit
            neigh{ii}(end + 1, 1) = perm(j); %#ok<AGROW>
            j = j + 1;
        end
    end
end

neigh = cellfun(@(v) uint32(v(:)), neigh, 'UniformOutput', false);
end

function [vPlus, vMinus, magPlus, magMinus] = local_gather_candidate_vectors(me, rpo, pair_me_global, pair_rpo_local)
%LOCAL_GATHER_CANDIDATE_VECTORS Gather the pair-level Moon-relative velocity vectors used by the flyby filter.

vPlus = me.vInf(double(pair_me_global), :);
vMinus = rpo.vInf(double(pair_rpo_local), :);
magPlus = me.vInfMag(double(pair_me_global));
magMinus = rpo.vInfMag(double(pair_rpo_local));
end

function ok = local_deflection_feasible(vPlus, vMinus, magPlus, magMinus, rpmin, muM)
%LOCAL_DEFLECTION_FEASIBLE Apply the Moon flyby deflection-feasibility filter.

dp = sum(vMinus .* vPlus, 2);
den = max(magMinus .* magPlus, realmin);
cosDelta = max(-1, min(1, dp ./ den));

s = 1 ./ (1 + (rpmin .* (magMinus .^ 2)) / muM);
s = max(0, min(1, s));
cosDeltaMax = 1 - 2 * (s .^ 2);
ok = (cosDelta >= cosDeltaMax);
end

function [chosenME_global, RPOMoonBest_idx, minDeltavp, batch_summary, filter_counts] = local_select_best_pairs_by_dvp(me, rpo, feasibleLimit, rpmin, muM, nSol, me_phase_batches, rpo_phase_batches, batch_summary)
%LOCAL_SELECT_BEST_PAIRS_BY_DVP Scan Moon-Earth rows in increasing Earth-flyby delta-v and stop once the best nSol are feasible.

chosenME_global = zeros(0, 1, 'uint32');
RPOMoonBest_idx = zeros(0, 1, 'uint32');
minDeltavp = zeros(0, 1);
pair_me_before_global = zeros(0, 1, 'uint32');
pair_rpo_before_local = zeros(0, 1, 'uint32');
pair_me_after_global = zeros(0, 1, 'uint32');
pair_rpo_after_local = zeros(0, 1, 'uint32');
me_candidate_global = vertcat(me_phase_batches{:});
if isempty(me_candidate_global)
    filter_counts = local_build_filter_counts(pair_me_before_global, pair_rpo_before_local, pair_me_after_global, pair_rpo_after_local);
    return
end

[~, order] = sort(me.dvpMag(double(me_candidate_global)), 'ascend');
me_candidate_global = me_candidate_global(order);

for k = 1:numel(me_candidate_global)
    me_global = me_candidate_global(k);
    batch_id = double(me.phase_sample_idx(double(me_global)));
    rpo_idx = rpo_phase_batches{batch_id};
    batch_summary.moonEarthCount(batch_id) = numel(me_phase_batches{batch_id});
    batch_summary.rpoCount(batch_id) = numel(rpo_idx);

    if isempty(rpo_idx)
        continue
    end

    % Search only inside the matched phase batch for close |vInf| neighbors.
    neigh = local_find_neighbors_by_magnitude(rpo.vInfMag(rpo_idx), me.vInfMag(double(me_global)), feasibleLimit);
    if isempty(neigh{1})
        continue
    end

    rpo_local = uint32(rpo_idx(double(neigh{1})));
    me_global_vec = repmat(uint32(me_global), numel(rpo_local), 1);
    pair_me_before_global = [pair_me_before_global; me_global_vec(:)]; %#ok<AGROW>
    pair_rpo_before_local = [pair_rpo_before_local; rpo_local(:)]; %#ok<AGROW>

    % Apply the same flyby-feasibility filter before accepting the row.
    [vPlus, vMinus, magPlus, magMinus] = local_gather_candidate_vectors(me, rpo, me_global_vec, rpo_local);
    ok = local_deflection_feasible(vPlus, vMinus, magPlus, magMinus, rpmin, muM);
    rpo_local = rpo_local(ok);
    batch_summary.pairCountBeforeFlyby(batch_id) = batch_summary.pairCountBeforeFlyby(batch_id) + numel(neigh{1});
    batch_summary.moonEarthRowsBeforeFlyby(batch_id) = batch_summary.moonEarthRowsBeforeFlyby(batch_id) + 1;
    batch_summary.rpoRowsBeforeFlyby(batch_id) = batch_summary.rpoRowsBeforeFlyby(batch_id) + numel(unique(double(rpo_idx(double(neigh{1})))));
    batch_summary.pairCountAfterFlyby(batch_id) = batch_summary.pairCountAfterFlyby(batch_id) + numel(rpo_local);
    batch_summary.pairCount(batch_id) = batch_summary.pairCountAfterFlyby(batch_id);
    if ~isempty(rpo_local)
        pair_me_after_global = [pair_me_after_global; repmat(uint32(me_global), numel(rpo_local), 1)]; %#ok<AGROW>
        pair_rpo_after_local = [pair_rpo_after_local; rpo_local(:)]; %#ok<AGROW>
        batch_summary.moonEarthRowsAfterFlyby(batch_id) = batch_summary.moonEarthRowsAfterFlyby(batch_id) + 1;
        batch_summary.rpoRowsAfterFlyby(batch_id) = batch_summary.rpoRowsAfterFlyby(batch_id) + numel(unique(double(rpo_local)));
    end

    if isempty(rpo_local)
        continue
    end

    % Keep the closest |vInf| neighbor for the selected Moon-Earth row.
    diff_mag = abs(rpo.vInfMag(double(rpo_local)) - me.vInfMag(double(me_global)));
    [~, best_idx] = min(diff_mag);
    chosenME_global(end + 1, 1) = uint32(me_global); %#ok<AGROW>
    RPOMoonBest_idx(end + 1, 1) = uint32(rpo_local(best_idx)); %#ok<AGROW>
    minDeltavp(end + 1, 1) = me.dvpMag(double(me_global)); %#ok<AGROW>

    if numel(chosenME_global) >= nSol
        break
    end
end

filter_counts = local_build_filter_counts(pair_me_before_global, pair_rpo_before_local, pair_me_after_global, pair_rpo_after_local);
end

function solutionTable = local_pack_solution_table(me, rpo, chosenME_global, RPOMoonBest_idx, minDeltavp, RPOMoonDatabase, timeSelection, geometryPolicy, phase_samples)
%LOCAL_PACK_SOLUTION_TABLE Pack the selected geometry-aware patch rows into the legacy output table plus extra metadata.

if isempty(chosenME_global)
    solutionTable = local_empty_solution_table();
    return
end
if nargin < 7 || isempty(timeSelection)
    timeSelection = struct();
end

deltavp = me.dvpvec(double(chosenME_global), :);
vpEarthFlyby = me.Vp(double(chosenME_global), :);
rpEarthFlyby = me.R2(double(chosenME_global), :);
T_lambert = me.tof(double(chosenME_global));
trueAnomalyMoonRPO = rpo.theta_raw(double(RPOMoonBest_idx));
trueAnomalyMoonRPORaw = rpo.theta_raw(double(RPOMoonBest_idx));
trueAnomalyMoon = me.theta(double(chosenME_global));
RI = me.RI(double(chosenME_global), :);
VI = me.VI(double(chosenME_global), :);
deltaVRPO = rpo.dv(double(RPOMoonBest_idx));
dv0CR3BP = rpo.dv0(double(RPOMoonBest_idx), :);
T_RPOMoon = rpo.t_event(double(RPOMoonBest_idx));
vInfPlusMoon = me.vInf(double(chosenME_global), :);
vInfMinusMoon = rpo.vInf(double(RPOMoonBest_idx), :);
correctedMoonPhaseRPO = rpo.theta_corrected(double(RPOMoonBest_idx));
moonEarthPhaseSampleIndex = double(me.phase_sample_idx(double(chosenME_global)));
moonEarthPhaseSampleValue = phase_samples(moonEarthPhaseSampleIndex);
selectedGeometryEpoch = repmat(geometryPolicy.initial_geometry_time, numel(chosenME_global), 1);

if ismember("orbitPhase", RPOMoonDatabase.Properties.VariableNames)
    phase = rpo.phase(double(RPOMoonBest_idx));
else
    phase = zeros(numel(chosenME_global), 1);
end
rpoRowIndex = rpo.row_index(double(RPOMoonBest_idx));
phaseDistanceRequest = rpo.phase_distance(double(RPOMoonBest_idx));
seedTimeMode = repmat(string("optimization_only"), numel(chosenME_global), 1);
seedToFLowerBoundDays = NaN(numel(chosenME_global), 1);
seedToFExactDays = NaN(numel(chosenME_global), 1);
seedToFCapDays = NaN(numel(chosenME_global), 1);
seedToFKeep = false(numel(chosenME_global), 1);
seedFlybyInboundDays = NaN(numel(chosenME_global), 1);
if isfield(timeSelection, 'time_mode')
    seedTimeMode = repmat(string(timeSelection.time_mode), numel(chosenME_global), 1);
end
if isfield(timeSelection, 'selected_tof_lower_bound_days')
    seedToFLowerBoundDays = timeSelection.selected_tof_lower_bound_days(:);
end
if isfield(timeSelection, 'selected_tof_exact_days')
    seedToFExactDays = timeSelection.selected_tof_exact_days(:);
end
if isfield(timeSelection, 'max_tof_days')
    seedToFCapDays = repmat(timeSelection.max_tof_days, numel(chosenME_global), 1);
end
if isfield(timeSelection, 'selected_tof_keep')
    seedToFKeep = logical(timeSelection.selected_tof_keep(:));
end
if isfield(timeSelection, 'selected_flyby_inbound_days')
    seedFlybyInboundDays = timeSelection.selected_flyby_inbound_days(:);
end

solutionTable = table( ...
    minDeltavp, deltavp, vpEarthFlyby, rpEarthFlyby, T_lambert, ...
    trueAnomalyMoon, RI, VI, deltaVRPO, dv0CR3BP, T_RPOMoon, ...
    trueAnomalyMoonRPO, vInfMinusMoon, vInfPlusMoon, phase, rpoRowIndex, phaseDistanceRequest, ...
    seedTimeMode, seedToFLowerBoundDays, seedToFExactDays, seedToFCapDays, seedToFKeep, seedFlybyInboundDays, ...
    correctedMoonPhaseRPO, moonEarthPhaseSampleIndex, moonEarthPhaseSampleValue, selectedGeometryEpoch, trueAnomalyMoonRPORaw, ...
    'VariableNames', {'deltavpMag','deltavp','vpEarthFlyby','rpEarthFlyby','T_lambert', ...
                      'trueAnomalyMoonFlyby','RI','VI', ...
                      'deltaVRPO','dv0CR3BP','T_RPOMoon', ...
                      'trueAnomalyMoonRPO','vInfMinusMoon','vInfPlusMoon','phase', ...
                      'rpoRowIndex','phaseDistanceRequest', ...
                      'seedTimeMode','seedToFLowerBoundDays','seedToFExactDays','seedToFCapDays','seedToFKeep','seedFlybyInboundDays', ...
                      'correctedMoonPhaseRPO','moonEarthPhaseSampleIndex','moonEarthPhaseSampleValue','selectedGeometryEpoch','trueAnomalyMoonRPORaw'});

% Attach candidate epochs when the Moon-Earth generator provides them.
if ~isempty(me.moonEpochNd)
    solutionTable.moonEpochNd = me.moonEpochNd(double(chosenME_global));
    solutionTable.earthEpochNd = me.earthEpochNd(double(chosenME_global));
    solutionTable.epochNodeIndex = me.epochNodeIndex(double(chosenME_global));
    solutionTable.earthPeriSampleIndex = me.earthPeriSampleIndex(double(chosenME_global));
end
end

function diagnostics = local_build_empty_diagnostics(batch_summary, geometryPolicy, vInfMoon_within_constraints)
%LOCAL_BUILD_EMPTY_DIAGNOSTICS Build one empty diagnostics package that still exposes the geometry-aware batch summary.

diagnostics = struct();
diagnostics.time_summary = table();
diagnostics.time_policy = struct();
diagnostics.refinement_pool_table = table();
diagnostics.refinement_pool_summary = table();
diagnostics.geometry_policy = geometryPolicy;
diagnostics.phase_samples = zeros(0, 1);
diagnostics.phase_batch_summary = batch_summary;
diagnostics.vInfMoon_within_constraints = vInfMoon_within_constraints;
end

function tbl = local_empty_solution_table()
%LOCAL_EMPTY_SOLUTION_TABLE Build one empty solution table with the full patched-seed schema.

tbl = RPO2DS_empty_patch_solution_table();
end

function filterSummary = local_build_filter_summary(me_total, rpo_total, me_prefilter, me_with_neighbors, pair_count_before, pair_count_after, unique_me_before, unique_rpo_before, unique_me_after, unique_rpo_after)
%LOCAL_BUILD_FILTER_SUMMARY Build one compact filter-survival summary reused by logs and tests.

filterSummary = table( ...
    me_total, rpo_total, me_prefilter, me_with_neighbors, ...
    pair_count_before, pair_count_after, ...
    unique_me_before, unique_rpo_before, ...
    unique_me_after, unique_rpo_after, ...
    'VariableNames', {'moonEarthRowsTotal','rpoRowsTotal','moonEarthRowsAfterVInfPrefilter','moonEarthRowsWithNeighbors', ...
    'candidatePairsBeforeFlyby','candidatePairsAfterFlyby', ...
    'uniqueMoonEarthRowsBeforeFlyby','uniqueRPORowsBeforeFlyby', ...
    'uniqueMoonEarthRowsAfterFlyby','uniqueRPORowsAfterFlyby'});
end

function filterCounts = local_build_filter_counts(pair_me_before_global, pair_rpo_before_local, pair_me_after_global, pair_rpo_after_local)
%LOCAL_BUILD_FILTER_COUNTS Build one compact filter-count package from the actual candidate-pair arrays.

filterCounts = struct();
filterCounts.moonEarthRowsWithNeighbors = numel(unique(double(pair_me_before_global)));
filterCounts.candidatePairsBeforeFlyby = numel(pair_me_before_global);
filterCounts.candidatePairsAfterFlyby = numel(pair_me_after_global);
filterCounts.uniqueMoonEarthRowsBeforeFlyby = numel(unique(double(pair_me_before_global)));
filterCounts.uniqueRPORowsBeforeFlyby = numel(unique(double(pair_rpo_before_local)));
filterCounts.uniqueMoonEarthRowsAfterFlyby = numel(unique(double(pair_me_after_global)));
filterCounts.uniqueRPORowsAfterFlyby = numel(unique(double(pair_rpo_after_local)));
end
