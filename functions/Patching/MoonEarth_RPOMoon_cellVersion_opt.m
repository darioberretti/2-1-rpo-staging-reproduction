function [solutionTable, diagnostics] = MoonEarth_RPOMoon_cellVersion_opt(MoonEarthDatabase, RPOMoonDatabase, nSol, feasibleLimit, hMinMoonFlyby, timePolicy)
%MOONEARTH_RPOMOON_CELLVERSION_OPT Match Moon-Earth patches to the RPO-Moon database.
%
% PROTOTYPE:
%   [solutionTable, diagnostics] = MoonEarth_RPOMoon_cellVersion_opt( ...
%       MoonEarthDatabase, RPOMoonDatabase, nSol, feasibleLimit, ...
%       hMinMoonFlyby, timePolicy)
%
% DESCRIPTION:
%   Prefilter lunar v-infinity magnitudes, find neighboring speeds, and
%   enforce minimum lunar-flyby deflection geometry. Rank surviving pairs
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
%
% OUTPUT:
%   solutionTable     Selected patch candidates, or the standard empty table.
%   diagnostics       Prefilter/flyby counts and timing/selection summaries.

if nargin < 5 || isempty(hMinMoonFlyby)
    hMinMoonFlyby = 100; % [km]
end
if nargin < 4 || isempty(feasibleLimit)
    feasibleLimit = 1e-3; % [km/s]
end
if nargin < 6 || isempty(timePolicy)
    timePolicy = struct();
end
if ~isfield(timePolicy, 'mode') || isempty(timePolicy.mode)
    timePolicy.mode = 'optimization_only';
end
if ~isfield(timePolicy, 'max_tof_days') || isempty(timePolicy.max_tof_days)
    timePolicy.max_tof_days = inf;
end

% Constants
persistent MU_MOON R_MOON
if isempty(MU_MOON)
    MU_MOON = astroConstants(20);
    R_MOON  = astroConstants(30);
end
muM = MU_MOON;
rpmin = R_MOON + hMinMoonFlyby;

% Extract arrays once (tables add overhead in tight loops)
me = extract_moon_earth(MoonEarthDatabase);
rpo = extract_rpo_moon(RPOMoonDatabase);
me_total = numel(me.vInfMag);
rpo_total = numel(rpo.vInfMag);

% Prefilter by vInf magnitude band
[idxME_all, diagnostics] = prefilter_by_vinf(me.vInfMag, rpo.vInfMag, feasibleLimit);
if isempty(idxME_all)
    diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(idxME_all), 0, 0, 0, 0, 0, 0, 0);
    diagnostics.refinement_pool_table = table();
    diagnostics.refinement_pool_summary = table();
    solutionTable = empty_solution_table();
    return
end

% Neighbor search: for each ME |vInf| find RPO neighbors within feasibleLimit
me_mags = me.vInfMag(idxME_all);
neigh = find_neighbors_by_magnitude(rpo.vInfMag, me_mags, feasibleLimit);

% Keep only ME rows with at least one neighbor
counts = cellfun('length', neigh);
hasAny = counts > 0;
if ~any(hasAny)
    diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(idxME_all), 0, 0, 0, 0, 0, 0, 0);
    diagnostics.refinement_pool_table = table();
    diagnostics.refinement_pool_summary = table();
    solutionTable = empty_solution_table();
    return
end

idxME = idxME_all(hasAny);
me_mags = me_mags(hasAny);
neigh = neigh(hasAny);
counts = counts(hasAny);

% Flatten candidate pairs (I,J)
I = repelem(uint32((1:numel(idxME))'), counts);
J = vertcat(neigh{:});
pair_count_before = numel(I);
unique_me_before = numel(unique(idxME(double(I))));
unique_rpo_before = numel(unique(double(J)));

% Deflection feasibility filter (vectorized)
[vPlus, vMinus, magPlus, magMinus] = gather_candidate_vectors(me, rpo, idxME, I, J);
[ok] = deflection_feasible(vPlus, vMinus, magPlus, magMinus, rpmin, muM);
I = I(ok);
J = J(ok);
pair_count_after = numel(I);
unique_me_after = numel(unique(idxME(double(I))));
unique_rpo_after = numel(unique(double(J)));
diagnostics.filter_summary = local_build_filter_summary(me_total, rpo_total, numel(idxME_all), numel(idxME), ...
    pair_count_before, pair_count_after, unique_me_before, unique_rpo_before, unique_me_after, unique_rpo_after);

if isempty(I)
    diagnostics.refinement_pool_table = table();
    diagnostics.refinement_pool_summary = table();
    solutionTable = empty_solution_table();
    return
end

% Keep the original delta-v-driven pair selection unless the capped time mode is active.
if strcmpi(timePolicy.mode, 'max_time')
    timeSelection = RPO2DS_time_capped_seed(me, rpo, idxME, I, J, nSol, timePolicy);
    if isempty(timeSelection.selected_me_global)
        diagnostics.refinement_pool_table = timeSelection.refinement_pool_table;
        diagnostics.refinement_pool_summary = timeSelection.pool_summary;
        solutionTable = empty_solution_table();
        return
    end
    diagnostics.time_summary = timeSelection.summary_table;
    diagnostics.time_policy = rmfield(timeSelection, 'summary_table');
    diagnostics.refinement_pool_table = timeSelection.refinement_pool_table;
    diagnostics.refinement_pool_summary = timeSelection.pool_summary;
    solutionTable = pack_solution_table( ...
        me, rpo, timeSelection.selected_me_global, timeSelection.selected_rpo_local, ...
        timeSelection.selected_deltavp_mag, RPOMoonDatabase, diagnostics.time_policy);
else
    % Choose best ME rows by smallest Earth-flyby delta-v.
    [chosenME_global, chosenLocalME, minDeltavp] = select_best_me_rows(me.dvpMag, idxME, I, nSol);
    if isempty(chosenME_global)
        solutionTable = empty_solution_table();
        return
    end

    % For each chosen ME row, pick the RPO neighbor with closest |vInf|.
    RPOMoonBest_idx = choose_best_rpo_neighbor(rpo.vInfMag, me.vInfMag, idxME, I, J, chosenLocalME);

    % Pack outputs using the unchanged legacy ranking.
    diagnostics.time_summary = table();
    diagnostics.time_policy = struct('time_mode', char(lower(string(timePolicy.mode))), ...
        'selected_tof_lower_bound_days', NaN(size(minDeltavp(:))), ...
        'selected_tof_exact_days', NaN(size(minDeltavp(:))), ...
        'selected_tof_keep', false(size(minDeltavp(:))), ...
        'selected_flyby_inbound_days', NaN(size(minDeltavp(:))), ...
        'max_tof_days', timePolicy.max_tof_days);
    diagnostics.refinement_pool_table = table();
    diagnostics.refinement_pool_summary = table();
    solutionTable = pack_solution_table(me, rpo, chosenME_global, RPOMoonBest_idx, minDeltavp, RPOMoonDatabase, diagnostics.time_policy);
end

end
% Helpers
function me = extract_moon_earth(MoonEarthDatabase)
%EXTRACT_MOON_EARTH Read Moon-Earth patch vectors, times, phases, and costs into numeric arrays.

    me = struct();
    me.vInfMag = MoonEarthDatabase.vInfMoonMag;
    me.vInf    = MoonEarthDatabase.vInfMoon;
    me.dvpvec  = MoonEarthDatabase.deltavp;
    me.dvpMag  = MoonEarthDatabase.deltavpMag;
    me.Vp      = MoonEarthDatabase.VF;   % keep existing behavior
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

function rpo = extract_rpo_moon(RPOMoonDatabase)
%EXTRACT_RPO_MOON Read RPO-Moon encounter velocities, burns, times, and phase metadata.

    rpo = struct();
    rpo.vInfMag = RPOMoonDatabase.vInfMagMoon;
    rpo.vInf    = RPOMoonDatabase.vInfMoon_RTH;
    rpo.dv      = RPOMoonDatabase.deltaV;
    rpo.dv0     = RPOMoonDatabase.deltaV_vect;
    rpo.t_event = RPOMoonDatabase.t_event;
    rpo.theta   = deg2rad(RPOMoonDatabase.thetaMoon_deg);
    rpo.row_index = (1:height(RPOMoonDatabase)).';
    if ismember("rpoRowIndexOriginal", RPOMoonDatabase.Properties.VariableNames)
        rpo.row_index = RPOMoonDatabase.rpoRowIndexOriginal;
    end
    rpo.phase_distance = NaN(height(RPOMoonDatabase), 1);
    if ismember("phaseDistanceRequest", RPOMoonDatabase.Properties.VariableNames)
        rpo.phase_distance = RPOMoonDatabase.phaseDistanceRequest;
    end
    if ismember("orbitPhase", RPOMoonDatabase.Properties.VariableNames)
        rpo.phase = RPOMoonDatabase.orbitPhase;
    else
        rpo.phase = 0;
    end
end

function [idxME_all, diagnostics] = prefilter_by_vinf(me_vInfMag, rpo_vInfMag, feasibleLimit)
%PREFILTER_BY_VINF Keep Moon-Earth speeds within the RPO database speed band and tolerance.

    lo = min(rpo_vInfMag) - feasibleLimit;
    hi = max(rpo_vInfMag) + feasibleLimit;
    keepME = (me_vInfMag >= lo) & (me_vInfMag <= hi);
    idxME_all = find(keepME);
    diagnostics = struct('vInfMoon_within_constraints', numel(idxME_all));
end

function neigh = find_neighbors_by_magnitude(rpo_vInfMag, me_mags, feasibleLimit)
%FIND_NEIGHBORS_BY_MAGNITUDE Find RPO speed neighbors with rangesearch or the equivalent direct scan.

    if exist('createns','file') == 2 && exist('rangesearch','file') == 2
        Mdl = createns(rpo_vInfMag(:), 'NSMethod','kdtree');
        neigh = rangesearch(Mdl, me_mags(:), feasibleLimit);
    else
        [rpo_sorted, perm] = sort(rpo_vInfMag(:));
        neigh = cell(numel(me_mags),1);
        j0 = 1;
        for ii = 1:numel(me_mags)
            x = me_mags(ii);
            while j0 <= numel(rpo_sorted) && rpo_sorted(j0) < x - feasibleLimit
                j0 = j0 + 1;
            end
            j = j0;
            while j <= numel(rpo_sorted) && rpo_sorted(j) <= x + feasibleLimit
                neigh{ii}(end+1,1) = perm(j);
                j = j + 1;
            end
        end
    end

    % Normalize neighbor lists to column vectors and uint32
    neigh = cellfun(@(v) uint32(v(:)), neigh, 'UniformOutput', false);
end

function [vPlus, vMinus, magPlus, magMinus] = gather_candidate_vectors(me, rpo, idxME, I, J)
%GATHER_CANDIDATE_VECTORS Gather incoming/outgoing lunar velocities for flattened candidate pairs.

    vPlus   = me.vInf(idxME(I), :);
    vMinus  = rpo.vInf(J, :);
    magPlus = me.vInfMag(idxME(I));
    magMinus= rpo.vInfMag(J);
end

function ok = deflection_feasible(vPlus, vMinus, magPlus, magMinus, rpmin, muM)
%DEFLECTION_FEASIBLE Check the lunar turning angle against the minimum-periapsis limit.

    dp = sum(vMinus.*vPlus, 2);
    den = max(magMinus.*magPlus, realmin);
    cosDelta = max(-1, min(1, dp./den));

    s = 1 ./ (1 + (rpmin .* (magMinus.^2))/muM);
    s = max(0, min(1, s));
    cosDeltaMax = 1 - 2*(s.^2);

    ok = (cosDelta >= cosDeltaMax);
end

function [chosenME_global, chosenLocalME, minDeltavp] = select_best_me_rows(me_dvpMag, idxME, I, nSol)
%SELECT_BEST_ME_ROWS Rank distinct Moon-Earth rows by perigee cost and retain the requested count.

    uniqueLocalME = unique(I, 'stable');
    meIdxLocal_to_global = idxME(uniqueLocalME);
    dvpFeas = me_dvpMag(meIdxLocal_to_global);

    nPick = min(nSol, numel(dvpFeas));
    if nPick == 0
        chosenME_global = [];
        chosenLocalME = [];
        minDeltavp = [];
        return
    end

    [minDeltavp, order] = mink(dvpFeas, nPick);
    chosenLocalME = uniqueLocalME(order);
    chosenME_global = meIdxLocal_to_global(order);
end

function RPOMoonBest_idx = choose_best_rpo_neighbor(rpo_vInfMag, me_vInfMag, idxME, I, J, chosenLocalME)
%CHOOSE_BEST_RPO_NEIGHBOR Choose the lowest-departure-cost RPO neighbor for each Earth-side row.

    diffMag_all = abs(rpo_vInfMag(J) - me_vInfMag(idxME(I)));
    isChosen = ismember(I, chosenLocalME);
    I2 = I(isChosen);
    J2 = J(isChosen);
    d2 = diffMag_all(isChosen);

    nPick = numel(chosenLocalME);
    RPOMoonBest_idx = zeros(nPick,1,'uint32');
    for t = 1:nPick
        loc = chosenLocalME(t);
        msk = (I2 == loc);
        [~, kmin] = min(d2(msk));
        Js = J2(msk);
        RPOMoonBest_idx(t) = Js(kmin);
    end
end

function solutionTable = pack_solution_table(me, rpo, chosenME_global, RPOMoonBest_idx, minDeltavp, RPOMoonDatabase, timeSelection)
%PACK_SOLUTION_TABLE Combine selected Earth-side and RPO-side rows into the shared patch schema.

    if isempty(chosenME_global)
        solutionTable = empty_solution_table();
        return
    end
    if nargin < 7 || isempty(timeSelection)
        timeSelection = struct();
    end

    deltavp         = me.dvpvec(chosenME_global, :);
    vpEarthFlyby    = me.Vp(chosenME_global, :);
    rpEarthFlyby    = me.R2(chosenME_global, :);
    T_lambert       = me.tof(chosenME_global);
    trueAnomalyMoonRPO = rpo.theta(RPOMoonBest_idx);
    trueAnomalyMoon = me.theta(chosenME_global);
    RI              = me.RI(chosenME_global, :);
    VI              = me.VI(chosenME_global, :);
    deltaVRPO       = rpo.dv(RPOMoonBest_idx);
    dv0CR3BP        = rpo.dv0(RPOMoonBest_idx, :);
    T_RPOMoon       = rpo.t_event(RPOMoonBest_idx);
    vInfPlusMoon    = me.vInf(chosenME_global, :);
    vInfMinusMoon   = rpo.vInf(RPOMoonBest_idx, :);

    if ismember("orbitPhase", RPOMoonDatabase.Properties.VariableNames)
        phase = rpo.phase(RPOMoonBest_idx);
    else
        phase = 0;
    end
    rpoRowIndex = rpo.row_index(RPOMoonBest_idx);
    phaseDistanceRequest = rpo.phase_distance(RPOMoonBest_idx);
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
        'VariableNames', {'deltavpMag','deltavp','vpEarthFlyby','rpEarthFlyby','T_lambert', ...
                          'trueAnomalyMoonFlyby','RI','VI', ...
                          'deltaVRPO','dv0CR3BP','T_RPOMoon', ...
                          'trueAnomalyMoonRPO','vInfMinusMoon','vInfPlusMoon','phase', ...
                          'rpoRowIndex','phaseDistanceRequest', ...
                          'seedTimeMode','seedToFLowerBoundDays','seedToFExactDays','seedToFCapDays','seedToFKeep','seedFlybyInboundDays'});

    % Attach candidate epochs when the Moon-Earth generator provides them.
    if ~isempty(me.moonEpochNd)
        solutionTable.moonEpochNd = me.moonEpochNd(chosenME_global);
        solutionTable.earthEpochNd = me.earthEpochNd(chosenME_global);
        solutionTable.epochNodeIndex = me.epochNodeIndex(chosenME_global);
        solutionTable.earthPeriSampleIndex = me.earthPeriSampleIndex(chosenME_global);
    end
end

function tbl = empty_solution_table()
%EMPTY_SOLUTION_TABLE Return the canonical zero-row patched-solution table.

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
