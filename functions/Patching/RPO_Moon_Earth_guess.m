function [FullSolution, hMin, hBest, diagnostics] = RPO_Moon_Earth_guess(vInfPlusEarth, RPOMoonDatabase, params)
%RPO_MOON_EARTH_GUESS Build a patched RPO-Moon-Earth initial guess.
%
% PROTOTYPE:
%   [FullSolution, hMin, hBest, diagnostics] = RPO_Moon_Earth_guess( ...
%       vInfPlusEarth, RPOMoonDatabase, params)
%
% DESCRIPTION:
%   Increase Earth-perigee altitude until the patch matcher finds a feasible
%   candidate, then probe higher altitudes for a lower perigee impulse. The
%   active time/geometry policy controls screening and ranking. This builds
%   an optimization seed; it does not solve the multiple-shooting NLP.
%
% INPUT:
%   vInfPlusEarth   [3 elements] Earth-relative outgoing ECI v-infinity [km/s].
%   RPOMoonDatabase  Historical or phase-aware RPO-Moon candidate table.
%   params           Struct with N, hStart [km], sweep_rad [rad], nSol,
%                      feasibleLimit [km/s], and optional nMoonOrbitSamples
%                      (default 100), time_mode/maxTOF [days], geometry_mode/
%                      initial_geometry_time [TU]. See RPO2DS_transfer.
%
% OUTPUT:
%   FullSolution     Selected patch table with source-row, velocity, flyby,
%                      phase, and timing fields consumed by RPO2DS_transfer.
%   hMin, hBest      First feasible and selected Earth-perigee altitudes [km].
%   diagnostics      Altitude scan, candidate screening, and timing metadata;
%                      an exhausted search returns an empty patch selection.

N = params.N;
hMinEarthFlyby = params.hStart;
minimum_height = hMinEarthFlyby;
sweep_rad = params.sweep_rad;
nSol = params.nSol;
feasibleLimit = params.feasibleLimit;
if ~isfield(params, 'nMoonOrbitSamples') || isempty(params.nMoonOrbitSamples)
    params.nMoonOrbitSamples = 100;
end
timePolicy = local_build_time_policy(params, nSol);
geometryPolicy = local_build_geometry_policy(params);

maxIter = 80;
iter = 0;
diagnostics = struct();
diagnostics.altitude_scan = table();
diagnostics.best_patch_summary = table();
diagnostics.best_patch_diagnostics = struct();
diagnostics.time_mode = char(timePolicy.mode);
diagnostics.max_tof_days = timePolicy.max_tof_days;
diagnostics.nMoonOrbitSamples = params.nMoonOrbitSamples;
progress_logs = local_progress_enabled(params);
local_progress_log(progress_logs, 'PATCH', 'starting altitude search: hStart %.6f km, N %d, Moon samples %d, RPO rows %d', hMinEarthFlyby, N, params.nMoonOrbitSamples, height(RPOMoonDatabase));
% Step 1: increase hMin until a feasible solution exists
sol = [];
while true
    local_progress_log(progress_logs, 'PATCH', 'altitude search iteration %d: testing h %.6f km', iter + 1, hMinEarthFlyby);
    [FullSolution, sol, patch_diag] = evaluate_patch(vInfPlusEarth, RPOMoonDatabase, N, hMinEarthFlyby, sweep_rad, nSol, feasibleLimit, params.nMoonOrbitSamples, timePolicy, geometryPolicy, params);
    diagnostics.altitude_scan = local_append_altitude_scan(diagnostics.altitude_scan, hMinEarthFlyby, sol);
    if isempty(sol)
        local_progress_log(progress_logs, 'PATCH', 'no feasible seed at h %.6f km; increasing altitude', hMinEarthFlyby);
        hMinEarthFlyby = hMinEarthFlyby * 1.1;
        iter = iter + 1;
        if iter >= maxIter
            hMin = hMinEarthFlyby;
            hBest = hMinEarthFlyby;
            diagnostics.best_patch_summary = table();
            diagnostics.best_patch_diagnostics = struct();
            local_progress_log(progress_logs, 'PATCH', 'altitude search reached max iterations without a feasible seed');
            return
        end
    else
        diagnostics.best_patch_summary = patch_diag.time_summary;
        diagnostics.best_patch_diagnostics = patch_diag;
        local_progress_log(progress_logs, 'PATCH', 'first feasible seed found at h %.6f km; deltavp %.6f km/s', hMinEarthFlyby, sol.dvpMag);
        break
    end
end

hMin = hMinEarthFlyby;
% Step 2: probe higher altitudes to look for lower delta-vp
hMinEarthFlybyNew = hMinEarthFlyby * 1.1;
local_progress_log(progress_logs, 'PATCH', 'probing higher altitude h %.6f km for a better seed', hMinEarthFlybyNew);
[FullSolution2, sol2, patch_diag2] = evaluate_patch(vInfPlusEarth, RPOMoonDatabase, N, hMinEarthFlybyNew, sweep_rad, nSol, feasibleLimit, params.nMoonOrbitSamples, timePolicy, geometryPolicy, params);
diagnostics.altitude_scan = local_append_altitude_scan(diagnostics.altitude_scan, hMinEarthFlybyNew, sol2);

if isempty(sol2)
    best = false;
    local_progress_log(progress_logs, 'PATCH', 'higher-altitude probe at h %.6f km produced no seed', hMinEarthFlybyNew);
else
    if local_is_better_candidate(sol2, sol, timePolicy)
        best = true;
        local_progress_log(progress_logs, 'PATCH', 'higher-altitude probe improved seed: h %.6f km, deltavp %.6f km/s', hMinEarthFlybyNew, sol2.dvpMag);
        FullSolution = FullSolution2;
        sol = sol2;
        diagnostics.best_patch_summary = patch_diag2.time_summary;
        diagnostics.best_patch_diagnostics = patch_diag2;
        while true
            hMinEarthFlybyNew = hMinEarthFlybyNew * 1.1;
            local_progress_log(progress_logs, 'PATCH', 'continuing altitude improvement scan at h %.6f km', hMinEarthFlybyNew);
            [FullSolution2, sol2, patch_diag2] = evaluate_patch(vInfPlusEarth, RPOMoonDatabase, N, hMinEarthFlybyNew, sweep_rad, nSol, feasibleLimit, params.nMoonOrbitSamples, timePolicy, geometryPolicy, params);
            diagnostics.altitude_scan = local_append_altitude_scan(diagnostics.altitude_scan, hMinEarthFlybyNew, sol2);
            if isempty(sol2) || ~local_is_better_candidate(sol2, sol, timePolicy)
                local_progress_log(progress_logs, 'PATCH', 'altitude improvement scan stopped at h %.6f km', hMinEarthFlybyNew);
                break
            end
            FullSolution = FullSolution2;
            sol = sol2;
            diagnostics.best_patch_summary = patch_diag2.time_summary;
            diagnostics.best_patch_diagnostics = patch_diag2;
            local_progress_log(progress_logs, 'PATCH', 'altitude improvement accepted: h %.6f km, deltavp %.6f km/s', hMinEarthFlybyNew, sol2.dvpMag);
        end
        hMinEarthFlybyBest = hMinEarthFlybyNew / 1.1;
        if hMinEarthFlybyBest < minimum_height
            hMinEarthFlybyBest = hMinEarthFlybyBest * 1.1;
        end
    else
        best = false;
        local_progress_log(progress_logs, 'PATCH', 'higher-altitude probe did not improve seed');
    end
end

if ~best
    hMinEarthFlybyBest = hMinEarthFlyby;
end

hBest = hMinEarthFlybyBest;
local_progress_log(progress_logs, 'PATCH', 'altitude search complete: first feasible h %.6f km, selected h %.6f km', hMin, hBest);

end
% Local helper
function [FullSolution, selectionMetrics, diagnostics] = evaluate_patch(vInfPlusEarth, RPOMoonDatabase, N, hMinEarthFlyby, sweep_rad, nSol, feasibleLimit, nMoonOrbitSamples, timePolicy, geometryPolicy, params)
%EVALUATE_PATCH Build the Earth-side candidates and run the active time/geometry patch matcher.

    progress_logs = local_progress_enabled(params);
    % Build the Moon-Earth database with the active target-frame timing law.
    local_progress_log(progress_logs, 'PATCH', 'building Moon-Earth database at h %.6f km', hMinEarthFlyby);
    if local_use_time_aware_rth(params)
        [MoonEarthDatabase, sampling_mode] = local_build_rth_time_aware_database( ...
            vInfPlusEarth, RPOMoonDatabase, N, hMinEarthFlyby, nMoonOrbitSamples, params, timePolicy, geometryPolicy);
    elseif isfield(timePolicy, 'fixed_arrival_epoch_time') && ~isempty(timePolicy.fixed_arrival_epoch_time)
        n_fixed_arrival_tof_samples = max(nMoonOrbitSamples, 25);
        MoonEarthDatabase = interplanetary_MoonEarth_patch_fixed_arrival(vInfPlusEarth, N, hMinEarthFlyby, sweep_rad, 0, ...
            n_fixed_arrival_tof_samples, timePolicy.fixed_arrival_epoch_time, timePolicy.max_tof_days, timePolicy.TU_days);
        sampling_mode = "fixed_arrival_epoch";
    else
        perigee_sampling_type = local_perigee_sampling_type(params);
        MoonEarthDatabase = interplanetary_MoonEarth_patch_opt(vInfPlusEarth, N, ...
            hMinEarthFlyby, sweep_rad, perigee_sampling_type, nMoonOrbitSamples);
        if perigee_sampling_type == 1
            sampling_mode = "retrograde_arc";
        else
            sampling_mode = "full_circle";
        end
    end
    local_progress_log(progress_logs, 'PATCH', 'Moon-Earth database built: mode %s, rows %d, RPO rows %d', char(sampling_mode), height(MoonEarthDatabase), height(RPOMoonDatabase));
    local_progress_log(progress_logs, 'PATCH', 'pairing Moon-Earth and RPO-Moon candidates');
    if strcmpi(geometryPolicy.mode, 'prescribed_geometry')
        [FullSolution, diagnostics] = MoonEarth_RPOMoon_cellVersion_opt_geometry(MoonEarthDatabase, RPOMoonDatabase, nSol, feasibleLimit, [], timePolicy, geometryPolicy);
    else
        [FullSolution, diagnostics] = MoonEarth_RPOMoon_cellVersion_opt(MoonEarthDatabase, RPOMoonDatabase, nSol, feasibleLimit, [], timePolicy);
    end
    local_log_patch_diagnostics(progress_logs, diagnostics, FullSolution);
    diagnostics.moon_earth_generation = struct( ...
        'earth_peri_samples', N, ...
        'moon_orbit_samples', nMoonOrbitSamples, ...
        'sampling_mode', sampling_mode, ...
        'rows_generated', height(MoonEarthDatabase));
    selectionMetrics = local_extract_selection_metrics(FullSolution, timePolicy);
end

function sampling_type = local_perigee_sampling_type(params)
%LOCAL_PERIGEE_SAMPLING_TYPE Select the historical retrograde arc only for callers that request it.

sampling_type = 0;
if isfield(params, 'perigee_sampling_type') && ...
        ~isempty(params.perigee_sampling_type)
    sampling_type = params.perigee_sampling_type;
end
if ~(isscalar(sampling_type) && ismember(sampling_type, [0, 1]))
    error('RPO_Moon_Earth_guess:PerigeeSamplingType', ...
        'perigee_sampling_type must be 0 (full circle) or 1 (retrograde arc).');
end
end

function timePolicy = local_build_time_policy(params, nSol)
%LOCAL_BUILD_TIME_POLICY Build the initial-guess time-policy package passed to the patch matcher.

timePolicy = struct();
timePolicy.mode = 'optimization_only';
timePolicy.max_tof_days = inf;
if isfield(params, 'time_mode') && ~isempty(params.time_mode)
    timePolicy.mode = lower(string(params.time_mode));
end
if isfield(params, 'maxTOF') && ~isempty(params.maxTOF)
    timePolicy.max_tof_days = params.maxTOF;
end
if isfield(params, 'fixed_arrival_epoch_time') && ~isempty(params.fixed_arrival_epoch_time)
    timePolicy.fixed_arrival_epoch_time = params.fixed_arrival_epoch_time;
    timePolicy.TU_days = params.TU_days_operational;
    timePolicy.initial_epoch_time = params.initial_epoch_time;
    timePolicy.initial_rpo_phase = params.initial_rpo_phase;
    timePolicy.fixed_arrival_phase_tolerance = params.fixed_arrival_phase_tolerance;
    timePolicy.phase_origin = params.phase_origin;
    timePolicy.phase_period = params.phase_period;
    timePolicy.phase_step = params.phase_step;
end
end

function tf = local_use_time_aware_rth(params)
%LOCAL_USE_TIME_AWARE_RTH Check whether the seed target must be rebuilt at each Earth-perigee epoch.

tf = isfield(params, 'target_frame') && strcmpi(string(params.target_frame), "earth_sun_rth");
end

function [MoonEarthDatabase, sampling_mode] = local_build_rth_time_aware_database( ...
    vInfPlusEarth, RPOMoonDatabase, N, hMinEarthFlyby, nMoonOrbitSamples, params, timePolicy, geometryPolicy)
%LOCAL_BUILD_RTH_TIME_AWARE_DATABASE Build the RTH-aware Moon-Earth table for standard or fixed-arrival seeds.

% Reuse the fixed-arrival timing grid when the final epoch is prescribed.
if isfield(timePolicy, 'fixed_arrival_epoch_time') && ~isempty(timePolicy.fixed_arrival_epoch_time)
    final_epoch = timePolicy.fixed_arrival_epoch_time;
    vInfEpoch = local_target_vinf_at_epoch(vInfPlusEarth, params, final_epoch);
    n_fixed_arrival_tof_samples = max(nMoonOrbitSamples, 25);
    MoonEarthDatabase = interplanetary_MoonEarth_patch_fixed_arrival(vInfEpoch, N, hMinEarthFlyby, params.sweep_rad, 0, ...
        n_fixed_arrival_tof_samples, final_epoch, timePolicy.max_tof_days, timePolicy.TU_days);
    sampling_mode = "rth_time_aware_fixed_arrival";
    return
end

% Build the standard free-final-time epoch grid from the active geometry mode.
epoch_nodes = local_standard_rth_epoch_nodes(RPOMoonDatabase, nMoonOrbitSamples, geometryPolicy);
MoonEarthDatabase = local_interplanetary_moon_earth_patch_time_aware(vInfPlusEarth, epoch_nodes, N, hMinEarthFlyby, params);
if strcmpi(string(geometryPolicy.mode), "prescribed_geometry")
    sampling_mode = "rth_time_aware_prescribed_geometry";
else
    sampling_mode = "rth_time_aware_full_circle";
end
end

function epoch_nodes = local_standard_rth_epoch_nodes(RPOMoonDatabase, nMoonOrbitSamples, geometryPolicy)
%LOCAL_STANDARD_RTH_EPOCH_NODES Build Moon-arrival epochs for the standard RTH seed search.

% Prescribed geometry knows the departure epoch, so each RPO row implies one Moon epoch.
if strcmpi(string(geometryPolicy.mode), "prescribed_geometry")
    epoch_nodes = geometryPolicy.initial_geometry_time + RPOMoonDatabase.t_event(:);
    epoch_nodes = unique(epoch_nodes, 'stable');
    return
end

% Free geometry keeps the legacy full-circle Moon sampling expressed as epochs.
epoch_nodes = linspace(0, 2 * pi, nMoonOrbitSamples + 1).';
epoch_nodes(end) = [];
end

function MoonEarthDatabase = local_interplanetary_moon_earth_patch_time_aware(vInfPlusEarth, epochNodesNd, N, hMinEarthFlyby, params)
%LOCAL_INTERPLANETARY_MOON_EARTH_PATCH_TIME_AWARE Build Moon-Earth rows whose outgoing Earth vInf target is epoch-aware.

% Store the patched-conic constants used by the Earth-Moon seed model.
rEM = 384400;
muE = 3.98600433e+05;
RE = astroConstants(23);
rp = RE + hMinEarthFlyby;
TU = 1 / 2.665e-6;
TU_days = TU / 86400;

% Preallocate the legacy columns plus explicit epoch metadata.
nEpoch = numel(epochNodesNd);
nRows = nEpoch * N;
progress_logs = local_progress_enabled(params);
local_progress_log(progress_logs, 'PATCH', 'RTH time-aware Moon-Earth build: epoch nodes %d, N %d, rows %d', nEpoch, N, nRows);
R1_all = zeros(nRows, 3);
R2_all = zeros(nRows, 3);
VI_all = zeros(nRows, 3);
VF_all = zeros(nRows, 3);
ToF_all = zeros(nRows, 1);
ToF_days_all = zeros(nRows, 1);
vInfMoon_all = zeros(nRows, 3);
vInfMoonMag_all = zeros(nRows, 1);
a_all = zeros(nRows, 1);
e_all = zeros(nRows, 1);
eps_all = zeros(nRows, 1);
theta_all = zeros(nRows, 1);
dvp_all = zeros(nRows, 3);
dvpMag_all = zeros(nRows, 1);
Vp_all = zeros(nRows, 3);
moonEpochNd = zeros(nRows, 1);
earthEpochNd = zeros(nRows, 1);
epochNodeIndex = zeros(nRows, 1);
earthPeriSampleIndex = zeros(nRows, 1);

% Sweep Moon epochs and Earth-perigee azimuth samples.
idx = 0;
progress_stride = max(1, ceil(nEpoch / 10));
for i = 1:nEpoch
    if i == 1 || mod(i, progress_stride) == 0 || i == nEpoch
        local_progress_log(progress_logs, 'PATCH', 'RTH time-aware Moon-Earth build progress: epoch %d/%d', i, nEpoch);
    end
    tMoon = epochNodesNd(i);
    thetaMoon = wrapTo2Pi(tMoon);
    R1 = rEM * [cos(thetaMoon), sin(thetaMoon), 0];
    for j = 1:N
        idx = idx + 1;
        az0 = 2 * pi * (j - 1) / N;
        tEarthGuess = tMoon + local_rth_seed_nominal_lambert_days(params) / TU_days;
        tEarthTarget = tEarthGuess;
        for it = 1:local_rth_seed_target_iterations(params)
            tEarthTarget = tEarthGuess;
            [A, E, VI, VF, ToF, R2, Vp] = local_time_aware_row( ...
                vInfPlusEarth, R1, rp, az0, tMoon, tEarthTarget, params, muE, TU);
            tEarthGuess = tMoon + ToF / TU;
        end

        % Rebuild the final row at the epoch that will be stored with the seed.
        tEarthTarget = tEarthGuess;
        [A, E, VI, VF, ToF, R2, Vp] = local_time_aware_row( ...
            vInfPlusEarth, R1, rp, az0, tMoon, tEarthTarget, params, muE, TU);

        % Convert Moon-departure velocity into the local Moon RTH frame.
        [vInfMoon, vInfMoonMag] = local_moon_vinf_from_vi(R1, VI, muE, rEM);
        deltavp = -VF(:).' + Vp(:).';

        % Store the row in the legacy Moon-Earth schema.
        R1_all(idx, :) = R1;
        R2_all(idx, :) = R2;
        VI_all(idx, :) = VI;
        VF_all(idx, :) = VF;
        ToF_all(idx) = ToF;
        ToF_days_all(idx) = ToF / 86400;
        vInfMoon_all(idx, :) = vInfMoon;
        vInfMoonMag_all(idx) = vInfMoonMag;
        a_all(idx) = A;
        e_all(idx) = E;
        eps_all(idx) = 0.5 * dot(VI, VI) - muE / rEM;
        theta_all(idx) = thetaMoon;
        dvp_all(idx, :) = deltavp;
        dvpMag_all(idx) = norm(deltavp);
        Vp_all(idx, :) = Vp;
        moonEpochNd(idx) = tMoon;
        earthEpochNd(idx) = tEarthTarget;
        epochNodeIndex(idx) = i;
        earthPeriSampleIndex(idx) = j;
    end
end
local_progress_log(progress_logs, 'PATCH', 'RTH time-aware Moon-Earth build complete');

% Return the legacy table shape with extra epoch columns for diagnostics and pairing.
MoonEarthDatabase = table( ...
    R1_all, R2_all, VI_all, VF_all, ToF_all, ToF_days_all, ...
    vInfMoon_all, vInfMoonMag_all, a_all, e_all, eps_all, theta_all, ...
    dvp_all, dvpMag_all, Vp_all, moonEpochNd, earthEpochNd, epochNodeIndex, earthPeriSampleIndex, ...
    'VariableNames', {'R1','R2','VI','VF','ToF','ToF_days', ...
    'vInfMoon','vInfMoonMag','a','e','eps','thetaMoon', ...
    'deltavp','deltavpMag','Vp','moonEpochNd','earthEpochNd','epochNodeIndex','earthPeriSampleIndex'});
end

function [A, E, VI, VF, ToF, R2, Vp] = local_time_aware_row( ...
    vInfPlusEarth, R1, rp, az0, tMoon, tEarthGuess, params, muE, TU)
%LOCAL_TIME_AWARE_ROW Build one epoch-aware Moon-Earth Lambert row.

% Build the outgoing Earth hyperbola from the target at the guessed Earth epoch.
vInfEpoch = local_target_vinf_at_epoch(vInfPlusEarth, params, tEarthGuess);
[Rp, VpGrid] = rp_circle_from_vinfplus(vInfEpoch, rp, muE, 1, az0);
R2 = Rp(1, :);
Vp = VpGrid(1, :);

% Solve for the Lambert time that makes the Earth-perigee velocity perpendicular.
r2norm = norm(R2);
tofSeed = max(1, (tEarthGuess - tMoon) * TU);
tofScale = pi * sqrt((rp + norm(R1))^3 / muE);
maxTofSec = local_rth_seed_max_lambert_days(params) * 86400;
[A, E, VI, VF, ToF] = RPO2DS_solve_perp_lambert_tof(R1, R2, muE, tofSeed, tofScale, r2norm, 1e-3, 60, maxTofSec);
end

function vInfEpoch = local_target_vinf_at_epoch(vInfPlusEarth, params, tEarth)
%LOCAL_TARGET_VINF_AT_EPOCH Convert the configured seed target into ECI at one Earth-perigee epoch.

% Keep the ECI target untouched when this helper is reused outside RTH mode.
if ~local_use_time_aware_rth(params)
    vInfEpoch = vInfPlusEarth(:);
    return
end

% Use the seed Sun phase for free-theta runs because the NLP optimizes it later.
theta_sun0 = params.theta_sun0_value;
if strcmpi(string(params.theta_sun0_mode), "free")
    theta_sun0 = params.theta_sun0_seed;
end

% Rotate the Earth-Sun RTH target into ECI at the candidate final epoch.
theta_sun_final = theta_sun0 + params.ws_sun_em_nd * tEarth;
vInfEpoch = earth_sun_rth_to_eci(params.vinf_target_rth_km_s, theta_sun_final, tEarth);
end

function [vInfMoon, vInfMoonMag] = local_moon_vinf_from_vi(R1, VI, muE, rEM)
%LOCAL_MOON_VINF_FROM_VI Convert the Moon-departure inertial velocity into Moon-relative RTH vInf.

rdir = R1 / norm(R1);
hdir = [0 0 1];
tdir = cross(hdir, rdir);
VI_RTH = [dot(VI', rdir), dot(VI', tdir), dot(VI', hdir)];
vMoon = [0 sqrt(muE / rEM) 0];
vInfMoon = VI_RTH(:).' - vMoon;
vInfMoonMag = norm(vInfMoon);
end

function nIter = local_rth_seed_target_iterations(params)
%LOCAL_RTH_SEED_TARGET_ITERATIONS Read the number of fixed-point updates used by the RTH seed builder.

nIter = 2;
if isfield(params, 'rth_seed_target_iterations') && ~isempty(params.rth_seed_target_iterations)
    nIter = params.rth_seed_target_iterations;
end
end

function tofDays = local_rth_seed_nominal_lambert_days(params)
%LOCAL_RTH_SEED_NOMINAL_LAMBERT_DAYS Read the nominal Moon-Earth Lambert duration seed in days.

tofDays = 8;
if isfield(params, 'rth_seed_nominal_lambert_tof_days') && ~isempty(params.rth_seed_nominal_lambert_tof_days)
    tofDays = params.rth_seed_nominal_lambert_tof_days;
end
end

function tofDays = local_rth_seed_max_lambert_days(params)
%LOCAL_RTH_SEED_MAX_LAMBERT_DAYS Read the maximum Moon-Earth Lambert duration used during seed construction.

tofDays = params.maxTOF;
if isfield(params, 'maxTransferTOF') && ~isempty(params.maxTransferTOF)
    tofDays = params.maxTransferTOF;
end
if isfield(params, 'rth_seed_max_lambert_tof_days') && ~isempty(params.rth_seed_max_lambert_tof_days)
    tofDays = params.rth_seed_max_lambert_tof_days;
end
end

function geometryPolicy = local_build_geometry_policy(params)
%LOCAL_BUILD_GEOMETRY_POLICY Build the geometry-selection package passed to the patch matcher.

geometryPolicy = struct();
geometryPolicy.mode = 'free_geometry';
geometryPolicy.initial_geometry_time = 0;
if isfield(params, 'geometry_mode') && ~isempty(params.geometry_mode)
    geometryPolicy.mode = lower(string(params.geometry_mode));
end
if isfield(params, 'initial_geometry_time') && ~isempty(params.initial_geometry_time)
    geometryPolicy.initial_geometry_time = params.initial_geometry_time;
end
end

function selectionMetrics = local_extract_selection_metrics(FullSolution, timePolicy)
%LOCAL_EXTRACT_SELECTION_METRICS Read the compact selection metrics used to compare altitude-search seeds.

selectionMetrics = [];
if ~RPO2DS_has_patch_solution(FullSolution)
    return
end

selectionMetrics = struct();
selectionMetrics.dvpMag = FullSolution.deltavpMag(1);
selectionMetrics.tofLowerBoundDays = NaN;
selectionMetrics.tofExactDays = NaN;
selectionMetrics.tofKeep = false;
if strcmpi(timePolicy.mode, 'max_time')
    if ismember('seedToFLowerBoundDays', FullSolution.Properties.VariableNames)
        selectionMetrics.tofLowerBoundDays = FullSolution.seedToFLowerBoundDays(1);
    end
    if ismember('seedToFExactDays', FullSolution.Properties.VariableNames)
        selectionMetrics.tofExactDays = FullSolution.seedToFExactDays(1);
    end
    if ismember('seedToFKeep', FullSolution.Properties.VariableNames)
        selectionMetrics.tofKeep = logical(FullSolution.seedToFKeep(1));
    end
end
end

function isBetter = local_is_better_candidate(candidateNew, candidateCurrent, timePolicy)
%LOCAL_IS_BETTER_CANDIDATE Compare two altitude-search candidates using the active time policy.

if isempty(candidateNew)
    isBetter = false;
    return
end
if isempty(candidateCurrent)
    isBetter = true;
    return
end

% Keep the legacy delta-v comparison when the capped-TOF policy is disabled.
if ~strcmpi(timePolicy.mode, 'max_time')
    isBetter = candidateNew.dvpMag < candidateCurrent.dvpMag;
    return
end

% Prefer exact under-cap seeds over exact over-cap seeds.
if candidateNew.tofKeep ~= candidateCurrent.tofKeep
    isBetter = candidateNew.tofKeep;
    return
end

% Rank under-cap seeds by delta-v first and exact TOF second.
if candidateNew.tofKeep
    if candidateNew.dvpMag ~= candidateCurrent.dvpMag
        isBetter = candidateNew.dvpMag < candidateCurrent.dvpMag;
    else
        isBetter = candidateNew.tofExactDays < candidateCurrent.tofExactDays;
    end
    return
end

% Rank over-cap seeds by exact TOF first and delta-v second.
if candidateNew.tofExactDays ~= candidateCurrent.tofExactDays
    isBetter = candidateNew.tofExactDays < candidateCurrent.tofExactDays;
else
    isBetter = candidateNew.dvpMag < candidateCurrent.dvpMag;
end
end

function altitudeScan = local_append_altitude_scan(altitudeScan, altitudeKm, selectionMetrics)
%LOCAL_APPEND_ALTITUDE_SCAN Append one compact altitude-search row for logs and tests.

if isempty(altitudeScan)
    altitudeScan = table('Size', [0, 5], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'logical'}, ...
        'VariableNames', {'hMinEarthFlybyKm', 'deltavpMag', 'tofLowerBoundDays', 'tofExactDays', 'tofKeep'});
end

if isempty(selectionMetrics)
    newRow = table(altitudeKm, NaN, NaN, NaN, false, ...
        'VariableNames', altitudeScan.Properties.VariableNames);
else
    newRow = table(altitudeKm, selectionMetrics.dvpMag, selectionMetrics.tofLowerBoundDays, selectionMetrics.tofExactDays, selectionMetrics.tofKeep, ...
        'VariableNames', altitudeScan.Properties.VariableNames);
end
altitudeScan = [altitudeScan; newRow];
end

function local_log_patch_diagnostics(enabled, diagnostics, FullSolution)
%LOCAL_LOG_PATCH_DIAGNOSTICS Print compact pairing diagnostics after one patching attempt.

if isfield(diagnostics, 'filter_summary') && istable(diagnostics.filter_summary) && ~isempty(diagnostics.filter_summary)
    fs = diagnostics.filter_summary(1, :);
    local_progress_log(enabled, 'PATCH', 'pairing filter summary: ME rows %s, RPO rows %s, pairs before %s, pairs after %s', ...
        local_table_value(fs, 'moonEarthRowsTotal'), local_table_value(fs, 'rpoRowsTotal'), ...
        local_table_value(fs, 'candidatePairsBeforeFlyby'), local_table_value(fs, 'candidatePairsAfterFlyby'));
end
if isfield(diagnostics, 'refinement_pool_summary') && istable(diagnostics.refinement_pool_summary) && ~isempty(diagnostics.refinement_pool_summary)
    ps = diagnostics.refinement_pool_summary(1, :);
    local_progress_log(enabled, 'PATCH', 'refinement pool summary: feasible pairs %s, selected pool %s, selected count %s', ...
        local_table_value(ps, 'totalFeasiblePairs'), local_table_value(ps, 'selectedPoolRole'), local_table_value(ps, 'selectedCount'));
end
if RPO2DS_has_patch_solution(FullSolution)
    local_progress_log(enabled, 'PATCH', 'pairing selected seed: rows %d, deltavp %.6f km/s', height(FullSolution), FullSolution.deltavpMag(1));
else
    local_progress_log(enabled, 'PATCH', 'pairing selected no seed');
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
