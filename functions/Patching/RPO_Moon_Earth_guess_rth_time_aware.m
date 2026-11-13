function [FullSolution, hMin, hBest, diagnostics] = RPO_Moon_Earth_guess_rth_time_aware(RPOMoonDatabase, params, phase_template, constants)
%RPO_MOON_EARTH_GUESS_RTH_TIME_AWARE Build RTH-aware operational patch seeds at candidate epochs.
%
% PROTOTYPE:
%   [FullSolution, hMin, hBest, diagnostics] = ...
%       RPO_Moon_Earth_guess_rth_time_aware(RPOMoonDatabase, params, ...
%       phase_template, constants)
%
% DESCRIPTION:
%   Expand RPO-Moon rows over admissible waiting periods, assign Moon-arrival
%   epoch nodes, and build Moon-Earth legs with the RTH target evaluated at
%   perigee. Increase Earth-perigee altitude until a feasible pair is found.
%
% INPUT:
%   RPOMoonDatabase  Phase-aware table containing encounter states/velocities.
%   params           Operational parameters from RPO2DS_transfer, including
%                      initial phase/epoch [TU], wait/time caps, RTH velocity
%                      [km/s], Sun phase [rad], and sampling settings.
%   phase_template   Phase origin/period [TU] from RPO2DS_phase_policy.
%   constants        Transfer constants with DU [km], VU [km/s], and TU_days.
%
% OUTPUT:
%   FullSolution     Selected epoch-matched patch table.
%   hMin, hBest      Feasible/selected Earth-perigee altitudes [km].
%   diagnostics      Sampling, epoch-assignment, pairing, and filter summaries.

hMin = params.hStart;
hBest = hMin;
maxIter = 80;
progress_logs = local_progress_enabled(params);
local_progress_log(progress_logs, 'PATCH-RTH', 'starting operational RTH time-aware patch search: hStart %.6f km, N %d', hMin, params.N);

% Expand each RPO-Moon row over the admissible onboard waiting times.
local_progress_log(progress_logs, 'PATCH-RTH', 'expanding RPO database over admissible waits');
RPOExpandedBase = local_expand_rpo_database(RPOMoonDatabase, params, phase_template, constants);
local_progress_log(progress_logs, 'PATCH-RTH', 'expanded RPO rows: %d from base rows %d', height(RPOExpandedBase), height(RPOMoonDatabase));

% Increase the Earth flyby altitude until the epoch-aware pairing has candidates.
for iter = 1:maxIter
    rpEarth = constants.Re + hBest;
    local_progress_log(progress_logs, 'PATCH-RTH', 'altitude iteration %d/%d: h %.6f km', iter, maxIter, hBest);

    % Build Moon-Earth epoch nodes over the expanded Moon-arrival window.
    RPOExpanded = RPOExpandedBase;
    nEpochNodes = local_epoch_node_count(params);
    tNodeMin = min(RPOExpanded.moonEpochNd);
    tNodeMax = max(RPOExpanded.moonEpochNd);
    epochNodesNd = linspace(tNodeMin, tNodeMax, nEpochNodes).';
    RPOExpanded.closestEpochNode = local_closest_epoch_node(RPOExpanded.moonEpochNd, epochNodesNd);
    RPOExpanded.closestEpochNd = epochNodesNd(RPOExpanded.closestEpochNode);
    RPOExpanded.epochMismatchDays = abs(RPOExpanded.moonEpochNd - RPOExpanded.closestEpochNd) * constants.TU_days;
    local_progress_log(progress_logs, 'PATCH-RTH', 'epoch nodes ready: %d nodes, moon epoch window [%.6f, %.6f] TU', nEpochNodes, tNodeMin, tNodeMax);

    % Build the time-aware Moon-Earth table with the active RTH target at perigee.
    local_progress_log(progress_logs, 'PATCH-RTH', 'building Moon-Earth database at h %.6f km', hBest);
    MoonEarthDatabase = local_build_moon_earth_database(epochNodesNd, params, constants, rpEarth);
    local_progress_log(progress_logs, 'PATCH-RTH', 'Moon-Earth database built: rows %d', height(MoonEarthDatabase));

    % Pair only rows assigned to the same epoch node and rank the feasible pool.
    local_progress_log(progress_logs, 'PATCH-RTH', 'pairing epoch-batched Moon-Earth and RPO rows');
    [FullSolution, diagnostics] = local_pair_and_pack(RPOExpanded, MoonEarthDatabase, RPOMoonDatabase, params, constants);
    local_log_rth_diagnostics(progress_logs, diagnostics, FullSolution);
    if RPO2DS_has_patch_solution(FullSolution)
        hMin = hBest;
        local_progress_log(progress_logs, 'PATCH-RTH', 'found feasible operational RTH seed at h %.6f km', hBest);
        break
    end
    hBest = hBest * 1.1;
    local_progress_log(progress_logs, 'PATCH-RTH', 'no seed found; increasing h to %.6f km', hBest);
end

diagnostics.moon_earth_generation = struct( ...
    'earth_peri_samples', params.N, ...
    'epoch_nodes', nEpochNodes, ...
    'rows_generated', height(MoonEarthDatabase), ...
    'hMinEarthFlybyKm', hMin);
diagnostics.rth_time_aware = struct( ...
    'expanded_rpo_rows', height(RPOExpanded), ...
    'median_epoch_assignment_days', median(RPOExpanded.epochMismatchDays), ...
    'max_epoch_assignment_days', max(RPOExpanded.epochMismatchDays));
local_progress_log(progress_logs, 'PATCH-RTH', 'operational RTH time-aware patch search complete: hMin %.6f km, hBest %.6f km', hMin, hBest);
end

function RPOExpanded = local_expand_rpo_database(RPOMoonDatabase, params, phase_template, constants)
%LOCAL_EXPAND_RPO_DATABASE Copy each RPO-Moon row over each allowed waiting time from the onboard state.

nRows = height(RPOMoonDatabase);
phaseSeedAll = RPOMoonDatabase.orbitPhase(:);
tEventAll = RPOMoonDatabase.t_event(:);
vinfMoonAll = RPOMoonDatabase.vInfMoon_RTH;
vinfMagAll = RPOMoonDatabase.vInfMagMoon(:);
deltaVAll = RPOMoonDatabase.deltaV(:);
deltaVVecAll = RPOMoonDatabase.deltaV_vect;
thetaMoonAll = deg2rad(RPOMoonDatabase.thetaMoon_deg(:));
waitCapNd = local_initial_state_wait_cap_nd(params);

% Count valid wait replicas before allocating the expanded arrays.
nMinAll = max(ceil((params.initial_rpo_phase - phaseSeedAll) / phase_template.phase_period), 0);
wait0All = phaseSeedAll + nMinAll * phase_template.phase_period - params.initial_rpo_phase;
nMaxAll = nMinAll + floor((waitCapNd - wait0All) / phase_template.phase_period);
replicaCount = max(nMaxAll - nMinAll + 1, 0);
nExpanded = sum(replicaCount);

rpoExpandedRow = zeros(nExpanded, 1);
rpoRowIndex = zeros(nExpanded, 1);
waitReplica = zeros(nExpanded, 1);
phase = zeros(nExpanded, 1);
waitTimeDays = zeros(nExpanded, 1);
departureEpochNd = zeros(nExpanded, 1);
moonEpochNd = zeros(nExpanded, 1);
vInfMoonMinus = zeros(nExpanded, 3);
vInfMoonMinusMag = zeros(nExpanded, 1);
deltaVRPO = zeros(nExpanded, 1);
dv0CR3BP = zeros(nExpanded, 3);
trueAnomalyMoonRPO = zeros(nExpanded, 1);
T_RPOMoon = zeros(nExpanded, 1);

% Fill contiguous blocks for each RPO-Moon database row.
writeStart = 1;
for k = 1:nRows
    countK = replicaCount(k);
    if countK == 0
        continue
    end

    idxRange = writeStart:(writeStart + countK - 1);
    nWait = (nMinAll(k):nMaxAll(k)).';
    phaseUnwrapped = phaseSeedAll(k) + nWait * phase_template.phase_period;
    waitNd = phaseUnwrapped - params.initial_rpo_phase;
    departureEpoch = params.initial_epoch_time + waitNd;
    moonEpoch = departureEpoch + tEventAll(k);

    rpoExpandedRow(idxRange) = idxRange(:);
    rpoRowIndex(idxRange) = k;
    waitReplica(idxRange) = nWait;
    phase(idxRange) = phaseUnwrapped;
    waitTimeDays(idxRange) = waitNd * constants.TU_days;
    departureEpochNd(idxRange) = departureEpoch;
    moonEpochNd(idxRange) = moonEpoch;
    vInfMoonMinus(idxRange, :) = repmat(vinfMoonAll(k, :), countK, 1);
    vInfMoonMinusMag(idxRange) = vinfMagAll(k);
    deltaVRPO(idxRange) = deltaVAll(k);
    dv0CR3BP(idxRange, :) = repmat(deltaVVecAll(k, :), countK, 1);
    trueAnomalyMoonRPO(idxRange) = thetaMoonAll(k);
    T_RPOMoon(idxRange) = tEventAll(k);
    writeStart = writeStart + countK;
end

RPOExpanded = table(rpoExpandedRow, rpoRowIndex, waitReplica, phase, waitTimeDays, ...
    departureEpochNd, moonEpochNd, vInfMoonMinus, vInfMoonMinusMag, deltaVRPO, ...
    dv0CR3BP, trueAnomalyMoonRPO, T_RPOMoon);
end

function MoonEarthDatabase = local_build_moon_earth_database(epochNodesNd, params, constants, rpEarth)
%LOCAL_BUILD_MOON_EARTH_DATABASE Build one time-aware Moon-Earth table over all epoch nodes.

nNodes = numel(epochNodesNd);
nRows = nNodes * params.N;
progress_logs = local_progress_enabled(params);
local_progress_log(progress_logs, 'PATCH-RTH', 'Moon-Earth row build start: epoch nodes %d, N %d, rows %d', nNodes, params.N, nRows);
R1 = zeros(nRows, 3);
R2 = zeros(nRows, 3);
VI = zeros(nRows, 3);
VF = zeros(nRows, 3);
Vp = zeros(nRows, 3);
vInfMoon = zeros(nRows, 3);
vInfMoonMag = zeros(nRows, 1);
ToF = zeros(nRows, 1);
thetaMoon = zeros(nRows, 1);
moonEpochNd = zeros(nRows, 1);
earthEpochNd = zeros(nRows, 1);
deltavp = zeros(nRows, 3);
deltavpMag = zeros(nRows, 1);
epochNodeIndex = zeros(nRows, 1);
earthPeriSampleIndex = zeros(nRows, 1);

idx = 0;
progress_stride = max(1, ceil(nNodes / 10));
for i = 1:nNodes
    if i == 1 || mod(i, progress_stride) == 0 || i == nNodes
        local_progress_log(progress_logs, 'PATCH-RTH', 'Moon-Earth row build progress: epoch %d/%d', i, nNodes);
    end
    tMoon = epochNodesNd(i);
    theta = wrapTo2Pi(tMoon);
    rMoon = constants.DU * [cos(theta), sin(theta), 0];

    for j = 1:params.N
        idx = idx + 1;
        az0 = 2 * pi * (j - 1) / params.N;
        tEarthGuess = tMoon + local_perigee_root_seed_days(params) / constants.TU_days;
        tEarthTarget = tEarthGuess;

        for it = 1:local_target_iterations(params)
            tEarthTarget = tEarthGuess;
            [rpVec, vpVec, tofSec, viVec, vfVec] = local_one_moon_earth_row( ...
                rMoon, rpEarth, az0, tMoon, tEarthTarget, params, constants);
            tEarthGuess = tMoon + tofSec / constants.TU;
        end

        % Rebuild the final row at the epoch that will be stored with the seed.
        tEarthTarget = tEarthGuess;
        [rpVec, vpVec, tofSec, viVec, vfVec] = local_one_moon_earth_row( ...
            rMoon, rpEarth, az0, tMoon, tEarthTarget, params, constants);

        [vinfMoonVec, vinfMoonMag] = local_moon_vinf_from_vi(rMoon, viVec, constants);
        R1(idx, :) = rMoon;
        R2(idx, :) = rpVec;
        VI(idx, :) = viVec;
        VF(idx, :) = vfVec;
        Vp(idx, :) = vpVec;
        vInfMoon(idx, :) = vinfMoonVec;
        vInfMoonMag(idx) = vinfMoonMag;
        ToF(idx) = tofSec;
        thetaMoon(idx) = theta;
        moonEpochNd(idx) = tMoon;
        earthEpochNd(idx) = tEarthTarget;
        deltavp(idx, :) = -vfVec + vpVec;
        deltavpMag(idx) = norm(deltavp(idx, :));
        epochNodeIndex(idx) = i;
        earthPeriSampleIndex(idx) = j;
    end
end
local_progress_log(progress_logs, 'PATCH-RTH', 'Moon-Earth row build complete');

MoonEarthDatabase = table(R1, R2, VI, VF, Vp, vInfMoon, vInfMoonMag, ToF, ...
    thetaMoon, moonEpochNd, earthEpochNd, deltavp, deltavpMag, epochNodeIndex, earthPeriSampleIndex);
end

function [rpVec, vpVec, tofSec, viVec, vfVec] = local_one_moon_earth_row(rMoon, rpEarth, az0, tMoon, tEarthGuess, params, constants)
%LOCAL_ONE_MOON_EARTH_ROW Build one Moon-Earth row for one epoch node and one Earth-perigee azimuth.

thetaSunFinal = params.theta_sun0_value + params.ws_sun_em_nd * tEarthGuess;
if strcmpi(string(params.theta_sun0_mode), "free")
    thetaSunFinal = params.theta_sun0_seed + params.ws_sun_em_nd * tEarthGuess;
end
vinfPlusEarth = earth_sun_rth_to_eci(params.vinf_target_rth_km_s, thetaSunFinal, tEarthGuess).';
[Rp, Vp] = rp_circle_from_vinfplus(vinfPlusEarth, rpEarth, constants.muE, 1, az0);
rpVec = Rp(1, :);
vpVec = Vp(1, :);
r2norm = norm(rpVec);
tofSeed = max(1, (tEarthGuess - tMoon) * constants.TU);
tofScale = pi * sqrt((rpEarth + constants.DU)^3 / constants.muE);
maxTofSec = local_perigee_root_max_days(params) * 86400;
[~, ~, viVec, vfVec, tofSec] = RPO2DS_solve_perp_lambert_tof(rMoon, rpVec, constants.muE, tofSeed, tofScale, r2norm, 1e-3, 60, maxTofSec);
viVec = viVec(:).';
vfVec = vfVec(:).';
end

function [FullSolution, diagnostics] = local_pair_and_pack(RPOExpanded, MoonEarthDatabase, RPOMoonDatabase, params, ~)
%LOCAL_PAIR_AND_PACK Pair Moon-Earth rows to expanded RPO-Moon rows inside each epoch batch.

rpmin = astroConstants(30) + 100;
nNodes = max(MoonEarthDatabase.epochNodeIndex);
progress_logs = local_progress_enabled(params);
local_progress_log(progress_logs, 'PATCH-RTH', 'epoch-batched pairing start: nodes %d, ME rows %d, expanded RPO rows %d', nNodes, height(MoonEarthDatabase), height(RPOExpanded));
pairRpo = cell(nNodes, 1);
pairMe = cell(nNodes, 1);
progress_stride = max(1, ceil(nNodes / 10));

for nodeIdx = 1:nNodes
    if nodeIdx == 1 || mod(nodeIdx, progress_stride) == 0 || nodeIdx == nNodes
        local_progress_log(progress_logs, 'PATCH-RTH', 'epoch-batched pairing progress: node %d/%d', nodeIdx, nNodes);
    end
    rpoIdx = find(RPOExpanded.closestEpochNode == nodeIdx);
    meIdx = find(MoonEarthDatabase.epochNodeIndex == nodeIdx);
    if isempty(rpoIdx) || isempty(meIdx)
        continue
    end

    neigh = local_find_neighbors_by_magnitude( ...
        RPOExpanded.vInfMoonMinusMag(rpoIdx), MoonEarthDatabase.vInfMoonMag(meIdx), params.feasibleLimit);
    counts = cellfun('length', neigh);
    hasAny = counts > 0;
    if ~any(hasAny)
        continue
    end

    localMe = find(hasAny);
    counts = counts(hasAny);
    neigh = neigh(hasAny);
    [neigh, counts] = local_trim_neighbors(neigh, counts, RPOExpanded.vInfMoonMinusMag(rpoIdx), ...
        MoonEarthDatabase.vInfMoonMag(meIdx), localMe, local_pair_cap(params));

    meLocalRep = repelem(localMe(:), counts);
    rpoLocalFlat = vertcat(neigh{:});
    pairRpo{nodeIdx} = rpoIdx(rpoLocalFlat);
    pairMe{nodeIdx} = meIdx(meLocalRep);
end

pairRpo = vertcat(pairRpo{:});
pairMe = vertcat(pairMe{:});
pair_count_before_empty = numel(pairRpo);
if isempty(pairRpo)
    FullSolution = table([], 'VariableNames', {'deltavpMag'});
    diagnostics = local_empty_diagnostics(height(MoonEarthDatabase), height(RPOExpanded), 0, 0);
    return
end

vOut = MoonEarthDatabase.vInfMoon(pairMe, :);
vIn = RPOExpanded.vInfMoonMinus(pairRpo, :);
magOut = MoonEarthDatabase.vInfMoonMag(pairMe);
magIn = RPOExpanded.vInfMoonMinusMag(pairRpo);
ok = local_deflection_feasible(vOut, vIn, magOut, magIn, rpmin, astroConstants(20));
pair_count_before_flyby = numel(ok);
pairRpo = pairRpo(ok);
pairMe = pairMe(ok);
if isempty(pairRpo)
    FullSolution = table([], 'VariableNames', {'deltavpMag'});
    diagnostics = local_empty_diagnostics(height(MoonEarthDatabase), height(RPOExpanded), pair_count_before_flyby, 0);
    return
end

moonVinfMismatch = abs(MoonEarthDatabase.vInfMoonMag(pairMe) - RPOExpanded.vInfMoonMinusMag(pairRpo));
rankTable = table(pairRpo, pairMe, MoonEarthDatabase.deltavpMag(pairMe), moonVinfMismatch, ...
    RPOExpanded.epochMismatchDays(pairRpo), ...
    'VariableNames', {'pairRpo','pairMe','deltavpMag','moonVinfMismatch','epochMismatchDays'});
rankTable = sortrows(rankTable, {'deltavpMag','moonVinfMismatch','epochMismatchDays'}, {'ascend','ascend','ascend'});
rankTable = rankTable(1:min(params.nSol, height(rankTable)), :);
FullSolution = local_pack_solution_table(rankTable, RPOExpanded, MoonEarthDatabase, RPOMoonDatabase);
diagnostics = local_empty_diagnostics(height(MoonEarthDatabase), height(RPOExpanded), pair_count_before_flyby, height(rankTable));
diagnostics.ranked_pairs = rankTable;
local_progress_log(progress_logs, 'PATCH-RTH', 'epoch-batched pairing complete: magnitude pairs %d, before flyby %d, after flyby %d', pair_count_before_empty, pair_count_before_flyby, height(rankTable));
end

function FullSolution = local_pack_solution_table(rankTable, RPOExpanded, MoonEarthDatabase, RPOMoonDatabase)
%LOCAL_PACK_SOLUTION_TABLE Pack the selected time-aware rows into the legacy FullSolution schema.

rpoIdx = rankTable.pairRpo;
meIdx = rankTable.pairMe;
minDeltavp = MoonEarthDatabase.deltavpMag(meIdx);
deltavp = MoonEarthDatabase.deltavp(meIdx, :);
vpEarthFlyby = MoonEarthDatabase.VF(meIdx, :);
rpEarthFlyby = MoonEarthDatabase.R2(meIdx, :);
T_lambert = MoonEarthDatabase.ToF(meIdx);
trueAnomalyMoonFlyby = MoonEarthDatabase.thetaMoon(meIdx);
RI = MoonEarthDatabase.R1(meIdx, :);
VI = MoonEarthDatabase.VI(meIdx, :);
deltaVRPO = RPOExpanded.deltaVRPO(rpoIdx);
dv0CR3BP = RPOExpanded.dv0CR3BP(rpoIdx, :);
T_RPOMoon = RPOExpanded.T_RPOMoon(rpoIdx);
trueAnomalyMoonRPO = RPOExpanded.trueAnomalyMoonRPO(rpoIdx);
vInfMinusMoon = RPOExpanded.vInfMoonMinus(rpoIdx, :);
vInfPlusMoon = MoonEarthDatabase.vInfMoon(meIdx, :);
phase = RPOExpanded.phase(rpoIdx);
rpoRowIndex = RPOExpanded.rpoRowIndex(rpoIdx);
phaseDistanceRequest = NaN(numel(rpoIdx), 1);
seedTimeMode = repmat("max_time", numel(rpoIdx), 1);
seedToFLowerBoundDays = NaN(numel(rpoIdx), 1);
seedToFExactDays = NaN(numel(rpoIdx), 1);
seedToFCapDays = NaN(numel(rpoIdx), 1);
seedToFKeep = false(numel(rpoIdx), 1);
seedFlybyInboundDays = NaN(numel(rpoIdx), 1);
correctedMoonPhaseRPO = wrapTo2Pi(RPOExpanded.departureEpochNd(rpoIdx) + T_RPOMoon);
moonEarthPhaseSampleIndex = MoonEarthDatabase.epochNodeIndex(meIdx);
moonEarthPhaseSampleValue = MoonEarthDatabase.moonEpochNd(meIdx);
selectedGeometryEpoch = RPOExpanded.departureEpochNd(rpoIdx);
trueAnomalyMoonRPORaw = trueAnomalyMoonRPO;
moonEpochNd = MoonEarthDatabase.moonEpochNd(meIdx);
earthEpochNd = MoonEarthDatabase.earthEpochNd(meIdx);
epochNodeIndex = MoonEarthDatabase.epochNodeIndex(meIdx);
earthPeriSampleIndex = MoonEarthDatabase.earthPeriSampleIndex(meIdx);

FullSolution = table( ...
    minDeltavp, deltavp, vpEarthFlyby, rpEarthFlyby, T_lambert, ...
    trueAnomalyMoonFlyby, RI, VI, deltaVRPO, dv0CR3BP, T_RPOMoon, ...
    trueAnomalyMoonRPO, vInfMinusMoon, vInfPlusMoon, phase, rpoRowIndex, phaseDistanceRequest, ...
    seedTimeMode, seedToFLowerBoundDays, seedToFExactDays, seedToFCapDays, seedToFKeep, seedFlybyInboundDays, ...
    correctedMoonPhaseRPO, moonEarthPhaseSampleIndex, moonEarthPhaseSampleValue, selectedGeometryEpoch, trueAnomalyMoonRPORaw, ...
    moonEpochNd, earthEpochNd, epochNodeIndex, earthPeriSampleIndex, ...
    'VariableNames', {'deltavpMag','deltavp','vpEarthFlyby','rpEarthFlyby','T_lambert', ...
    'trueAnomalyMoonFlyby','RI','VI', ...
    'deltaVRPO','dv0CR3BP','T_RPOMoon', ...
    'trueAnomalyMoonRPO','vInfMinusMoon','vInfPlusMoon','phase', ...
    'rpoRowIndex','phaseDistanceRequest', ...
    'seedTimeMode','seedToFLowerBoundDays','seedToFExactDays','seedToFCapDays','seedToFKeep','seedFlybyInboundDays', ...
    'correctedMoonPhaseRPO','moonEarthPhaseSampleIndex','moonEarthPhaseSampleValue','selectedGeometryEpoch','trueAnomalyMoonRPORaw', ...
    'moonEpochNd','earthEpochNd','epochNodeIndex','earthPeriSampleIndex'});

if ismember("rpoRowIndexOriginal", RPOMoonDatabase.Properties.VariableNames)
    FullSolution.rpoRowIndex = RPOMoonDatabase.rpoRowIndexOriginal(FullSolution.rpoRowIndex);
end
end

function diagnostics = local_empty_diagnostics(meRows, rpoRows, pairsBeforeFlyby, pairsAfterFlyby)
%LOCAL_EMPTY_DIAGNOSTICS Build compact diagnostics for the time-aware RTH patcher.

diagnostics = struct();
diagnostics.time_summary = table();
diagnostics.time_policy = struct('time_mode', 'max_time');
diagnostics.refinement_pool_table = table();
diagnostics.refinement_pool_summary = table();
diagnostics.filter_summary = table(meRows, rpoRows, pairsBeforeFlyby, pairsAfterFlyby, ...
    'VariableNames', {'moonEarthRowsTotal','rpoRowsTotal','candidatePairsBeforeFlyby','candidatePairsAfterFlyby'});
end

function neigh = local_find_neighbors_by_magnitude(rpoMag, meMag, feasibleLimit)
%LOCAL_FIND_NEIGHBORS_BY_MAGNITUDE Find local RPO magnitude neighbors for each Moon-Earth magnitude.

[rpoSorted, perm] = sort(rpoMag(:));
neigh = cell(numel(meMag), 1);
for k = 1:numel(meMag)
    lo = local_lower_bound(rpoSorted, meMag(k) - feasibleLimit);
    hi = local_upper_bound(rpoSorted, meMag(k) + feasibleLimit);
    if hi >= lo
        neigh{k} = perm(lo:hi);
    else
        neigh{k} = zeros(0, 1);
    end
end
end

function [neigh, counts] = local_trim_neighbors(neigh, counts, rpoMag, meMag, localMe, pairCap)
%LOCAL_TRIM_NEIGHBORS Keep the closest magnitude matches when one epoch node is too dense.

rawCount = sum(counts);
if rawCount <= pairCap
    return
end

keepEach = max(1, floor(pairCap / max(numel(neigh), 1)));
for k = 1:numel(neigh)
    idx = neigh{k};
    if numel(idx) <= keepEach
        continue
    end
    diffMag = abs(rpoMag(idx) - meMag(localMe(k)));
    [~, order] = mink(diffMag, keepEach);
    neigh{k} = idx(order);
    counts(k) = keepEach;
end
end

function ok = local_deflection_feasible(vInfOut, vInfIn, magOut, magIn, rpmin, muM)
%LOCAL_DEFLECTION_FEASIBLE Check vectorized ballistic Moon flyby rotation feasibility.

cosDelta = sum(vInfIn .* vInfOut, 2) ./ max(magIn .* magOut, realmin);
cosDelta = max(-1, min(1, cosDelta));
s = 1 ./ (1 + (rpmin .* magIn .^ 2) ./ muM);
cosDeltaMax = 1 - 2 .* s .^ 2;
ok = cosDelta >= cosDeltaMax;
end

function [vInfMoon, vInfMoonMag] = local_moon_vinf_from_vi(rMoon, viVec, constants)
%LOCAL_MOON_VINF_FROM_VI Convert inertial Moon-departure velocity into local Moon-relative RTH vInf.

rdir = rMoon / norm(rMoon);
hdir = [0 0 1];
tdir = cross(hdir, rdir);
VI_RTH = [dot(viVec, rdir), dot(viVec, tdir), dot(viVec, hdir)];
vMoon_RTH = [0 sqrt(constants.muE / constants.DU) 0];
vInfMoon = VI_RTH - vMoon_RTH;
vInfMoonMag = norm(vInfMoon);
end

function closestIdx = local_closest_epoch_node(tValues, epochNodesNd)
%LOCAL_CLOSEST_EPOCH_NODE Assign each expanded RPO-Moon epoch to the closest Moon-Earth epoch node.

closestIdx = zeros(numel(tValues), 1);
for k = 1:numel(tValues)
    [~, closestIdx(k)] = min(abs(epochNodesNd - tValues(k)));
end
end

function idx = local_lower_bound(sortedValues, value)
%LOCAL_LOWER_BOUND Return the first sorted index with value greater than or equal to target.

lo = 1;
hi = numel(sortedValues) + 1;
while lo < hi
    mid = floor((lo + hi) / 2);
    if sortedValues(mid) < value
        lo = mid + 1;
    else
        hi = mid;
    end
end
idx = lo;
end

function idx = local_upper_bound(sortedValues, value)
%LOCAL_UPPER_BOUND Return the last sorted index with value less than or equal to target.

lo = 1;
hi = numel(sortedValues) + 1;
while lo < hi
    mid = floor((lo + hi) / 2);
    if sortedValues(mid) <= value
        lo = mid + 1;
    else
        hi = mid;
    end
end
idx = lo - 1;
end

function waitCapNd = local_initial_state_wait_cap_nd(params)
%LOCAL_INITIAL_STATE_WAIT_CAP_ND Use the same two-RPO-period wait cap enforced by the NLP.

waitCapNd = 2 * (-pi / params.ws_sun_em_nd);
end

function nEpochNodes = local_epoch_node_count(params)
%LOCAL_EPOCH_NODE_COUNT Read the number of time-aware Moon-Earth epoch nodes.

nEpochNodes = params.nMoonOrbitSamples;
if isfield(params, 'rth_seed_epoch_nodes') && ~isempty(params.rth_seed_epoch_nodes)
    nEpochNodes = params.rth_seed_epoch_nodes;
end
end

function nIter = local_target_iterations(params)
%LOCAL_TARGET_ITERATIONS Read the fixed-point iterations for the active RTH target epoch.

nIter = 2;
if isfield(params, 'rth_seed_target_iterations') && ~isempty(params.rth_seed_target_iterations)
    nIter = params.rth_seed_target_iterations;
end
end

function tofDays = local_perigee_root_seed_days(params)
%LOCAL_PERIGEE_ROOT_SEED_DAYS Read the nominal Lambert seed used by the perigee time search.

tofDays = 8;
if isfield(params, 'rth_seed_nominal_lambert_tof_days') && ~isempty(params.rth_seed_nominal_lambert_tof_days)
    tofDays = params.rth_seed_nominal_lambert_tof_days;
end
end

function tofDays = local_perigee_root_max_days(params)
%LOCAL_PERIGEE_ROOT_MAX_DAYS Read the maximum Lambert TOF allowed during the perigee time search.

tofDays = params.maxTransferTOF;
if isfield(params, 'rth_seed_max_lambert_tof_days') && ~isempty(params.rth_seed_max_lambert_tof_days)
    tofDays = params.rth_seed_max_lambert_tof_days;
end
end

function pairCap = local_pair_cap(params)
%LOCAL_PAIR_CAP Read the maximum magnitude-filtered pair count stored per epoch node.

pairCap = max(5e4, 100 * params.nSol);
if isfield(params, 'rth_seed_max_pairs_per_epoch_node') && ~isempty(params.rth_seed_max_pairs_per_epoch_node)
    pairCap = params.rth_seed_max_pairs_per_epoch_node;
end
end

function local_log_rth_diagnostics(enabled, diagnostics, FullSolution)
%LOCAL_LOG_RTH_DIAGNOSTICS Print compact diagnostics for one operational RTH patching attempt.

if isfield(diagnostics, 'filter_summary') && istable(diagnostics.filter_summary) && ~isempty(diagnostics.filter_summary)
    fs = diagnostics.filter_summary(1, :);
    local_progress_log(enabled, 'PATCH-RTH', 'filter summary: ME rows %s, RPO rows %s, pairs before flyby %s, pairs after flyby %s', ...
        local_table_value(fs, 'moonEarthRowsTotal'), local_table_value(fs, 'rpoRowsTotal'), ...
        local_table_value(fs, 'candidatePairsBeforeFlyby'), local_table_value(fs, 'candidatePairsAfterFlyby'));
end
if RPO2DS_has_patch_solution(FullSolution)
    local_progress_log(enabled, 'PATCH-RTH', 'selected seed rows %d, deltavp %.6f km/s', height(FullSolution), FullSolution.deltavpMag(1));
else
    local_progress_log(enabled, 'PATCH-RTH', 'no selected seed in this attempt');
end
end

function value_text = local_table_value(row, variable_name)
%LOCAL_TABLE_VALUE Convert one optional scalar table value into a printable string.

if ~ismember(variable_name, row.Properties.VariableNames)
    value_text = 'n/a';
    return
end
value = row.(variable_name)(1);
if isnumeric(value) || islogical(value)
    value_text = num2str(value);
elseif isstring(value)
    value_text = char(value);
elseif iscell(value)
    value_text = char(string(value{1}));
else
    value_text = char(string(value));
end
end

function log_config = local_progress_enabled(params)
%LOCAL_PROGRESS_ENABLED Read the optional console progress-log flag from the params struct.

log_config = false;
if isfield(params, 'progress_logs') && ~isempty(params.progress_logs)
    log_config = params.progress_logs;
end
end

function local_progress_log(log_config, component, message, varargin)
%LOCAL_PROGRESS_LOG Print one timestamped progress message and flush the MATLAB console.

[enabled, t_start] = local_progress_state(log_config);
if ~enabled
    return
end
time_text = char(datetime('now', 'Format', 'HH:mm:ss'));
fprintf('[%s | +%.1fs] %s | %s\n', time_text, toc(t_start), component, sprintf(message, varargin{:}));
drawnow('limitrate');
end

function [enabled, t_start] = local_progress_state(log_config)
%LOCAL_PROGRESS_STATE Resolve boolean and elapsed-time origin for console progress logs.

enabled = false;
t_start = tic;
if isstruct(log_config)
    if isfield(log_config, 'enabled') && ~isempty(log_config.enabled)
        enabled = logical(log_config.enabled);
    end
    if isfield(log_config, 't_start') && ~isempty(log_config.t_start)
        t_start = log_config.t_start;
    end
else
    enabled = ~isempty(log_config) && logical(log_config);
end
if ~isstruct(log_config) && enabled && isappdata(0, 'RPO_progress_log_t_start')
    t_start = getappdata(0, 'RPO_progress_log_t_start');
end
end
