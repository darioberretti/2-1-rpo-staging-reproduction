function interplanetaryEarthMoon = interplanetary_MoonEarth_patch_fixed_arrival(vInfPlusEarth, N, hMinEarthFlyby, sweep_rad, type, nTofSamples, finalEpochTU, maxLambertTofDays, TU_days)
%INTERPLANETARY_MOONEARTH_PATCH_FIXED_ARRIVAL Build Moon-Earth candidates at a fixed perigee epoch.
%
% PROTOTYPE:
%   database = interplanetary_MoonEarth_patch_fixed_arrival(vInfPlusEarth, ...
%       N, hMinEarthFlyby, sweep_rad, type, nTofSamples, finalEpochTU, ...
%       maxLambertTofDays, TU_days)
%
% DESCRIPTION:
%   Fix the final Earth-perigee epoch and sample flight times backward to
%   recover the Moon departure geometry. Preserve the same patch table
%   columns used by the free-geometry construction.
%
% INPUT:
%   vInfPlusEarth      Outgoing Earth-relative ECI v-infinity [km/s].
%   N                  Number of Earth-perigee samples.
%   hMinEarthFlyby     Minimum Earth-perigee altitude [km].
%   sweep_rad, type    Perigee arc span [rad] and full/retrograde sampling flag.
%   nTofSamples        Number of Lambert flight-time samples.
%   finalEpochTU       Fixed Earth-perigee epoch [TU].
%   maxLambertTofDays   Maximum sampled flight duration [days].
%   TU_days            Days per normalized time unit.
%
% OUTPUT:
%   interplanetaryEarthMoon  Patch candidate table with positions [km],
%                             velocities [km/s], flight times [s/days],
%                             and final/departure epoch metadata.

if nargin < 6 || isempty(nTofSamples)
    nTofSamples = 100;
end

% Store the constants used by the patched-conic Moon-Earth model.
rEM = 384400;               % Earth-Moon distance [km]
muE = 3.98600433e+05;       % Earth gravitational parameter [km^3/s^2]
RE = astroConstants(23);    % Earth radius [km]
rp = RE + hMinEarthFlyby;   % Earth flyby radius [km]
day = 3600 * 24;            % Seconds per day [s]

% Keep the target outgoing infinity vector in row-vector form.
vinf_plus = vInfPlusEarth(:).';

% Build the sampled Earth-perigee circle exactly as in the legacy database.
if ~type
    [Rp, Vp] = rp_circle_from_vinfplus(vinf_plus, rp, muE, N);
else
    rM_vec = [rEM 0 0]';
    vM_vec = [0 sqrt(muE / rEM) 0]';
    arc = retrograde_rp_arc(vinf_plus, rp, muE, rM_vec, vM_vec, sweep_rad, N);
    if ~isempty(arc)
        Rp = arc.Rp;
        Vp = arc.Vp;
    else
        [Rp, Vp] = rp_circle_from_vinfplus(vinf_plus, rp, muE, N);
    end
end

% Sample Lambert durations inside the active fixed-arrival transfer window.
tof_min_days = max(1 / day, 1e-6);
tof_max_days = max(tof_min_days, maxLambertTofDays);
tof_days_grid = linspace(tof_min_days, tof_max_days, nTofSamples);
tof_sec_grid = tof_days_grid * day;

% Preallocate one flat table with N flyby samples times nTofSamples durations.
M = N * nTofSamples;
R1_all = zeros(M, 3);
R2_all = zeros(M, 3);
VI_all = zeros(M, 3);
VF_all = zeros(M, 3);
ToF_all = zeros(M, 1);
ToF_days_all = zeros(M, 1);
vInfMoon_all = zeros(M, 3);
vInfMoonMag_all = zeros(M, 1);
a_all = zeros(M, 1);
e_all = zeros(M, 1);
eps_all = zeros(M, 1);
theta_all = zeros(M, 1);
dvp_all = zeros(M, 3);
dvpMag_all = zeros(M, 1);
Vp_all = zeros(M, 3);
moonEpochNd_all = zeros(M, 1);
earthEpochNd_all = zeros(M, 1);
epochNodeIndex_all = zeros(M, 1);
earthPeriSampleIndex_all = zeros(M, 1);

% Try every Earth flyby sample against every fixed-arrival Lambert duration.
for j = 1:N
    R2 = Rp(j, :);
    vp = Vp(j, :);
    base = (j - 1) * nTofSamples;

    for k = 1:nTofSamples
        idx = base + k;
        ToF_star = tof_sec_grid(k);
        moon_epoch_tu = finalEpochTU - tof_days_grid(k) / TU_days;
        thetaMoon = wrapTo2Pi(moon_epoch_tu);
        R1 = rEM * [cos(thetaMoon), sin(thetaMoon), 0];

        % Solve the Moon-to-Earth Lambert leg for the implied Moon position.
        % Use the repository-owned FEX-Lambert wrapper; the old lambertMR
        % dependency is intentionally not redistributed.
        [VI, VF] = rpo2ds.solve_lambert(R1, R2, tof_days_grid(k), muE);
        [A, E] = local_orbit_elements(R1, VI, muE);
        VI = VI(:).';
        VF = VF(:).';

        % Convert the Moon-departure velocity into the local Moon RTH frame.
        eps_value = 0.5 * dot(VI, VI) - muE / rEM;
        rdir = R1 / norm(R1);
        hdir = [0 0 1]';
        tdir = cross(hdir, rdir);
        vMoon = [0 sqrt(muE / rEM) 0];
        VI_RTH = [dot(VI', rdir), dot(VI', tdir), dot(VI', hdir)];
        vInfMoon = VI_RTH(:).' - vMoon;
        vInfMoonMag = norm(vInfMoon);

        % Compute the perigee maneuver that turns the incoming leg into vInf+.
        deltavp = -VF(:).' + vp(:).';
        deltavpMag = norm(deltavp);

        % Store the candidate in the legacy Moon-Earth table schema.
        R1_all(idx, :) = R1;
        R2_all(idx, :) = R2;
        VI_all(idx, :) = VI;
        VF_all(idx, :) = VF;
        ToF_all(idx) = ToF_star;
        ToF_days_all(idx) = tof_days_grid(k);
        vInfMoon_all(idx, :) = vInfMoon;
        vInfMoonMag_all(idx) = vInfMoonMag;
        a_all(idx) = A;
        e_all(idx) = E;
        eps_all(idx) = eps_value;
        theta_all(idx) = thetaMoon;
        dvp_all(idx, :) = deltavp;
        dvpMag_all(idx) = deltavpMag;
        Vp_all(idx, :) = vp;
        moonEpochNd_all(idx) = moon_epoch_tu;
        earthEpochNd_all(idx) = finalEpochTU;
        epochNodeIndex_all(idx) = k;
        earthPeriSampleIndex_all(idx) = j;
    end
end

% Build the output table expected by the existing RPO-Moon patch matcher.
MoonEarthDatabase = table( ...
    R1_all, R2_all, VI_all, VF_all, ToF_all, ToF_days_all, ...
    vInfMoon_all, vInfMoonMag_all, ...
    a_all, e_all, eps_all, theta_all, dvp_all, dvpMag_all, Vp_all, ...
    moonEpochNd_all, earthEpochNd_all, epochNodeIndex_all, earthPeriSampleIndex_all, ...
    'VariableNames', {'R1','R2','VI','VF','ToF','ToF_days', ...
                      'vInfMoon','vInfMoonMag','a','e','eps','thetaMoon', ...
                      'deltavp','deltavpMag','Vp', ...
                      'moonEpochNd','earthEpochNd','epochNodeIndex','earthPeriSampleIndex'} );

% Remove numerical dust from vector-valued columns before returning.
zeroTol = 1e-10;
for v = intersect(["R1","R2","VI","VF","ToF","vInfMoon"], string(MoonEarthDatabase.Properties.VariableNames))
    MoonEarthDatabase.(v)(abs(MoonEarthDatabase.(v)) < zeroTol) = 0;
end

interplanetaryEarthMoon = MoonEarthDatabase;
end

function [semi_major_axis_km, eccentricity] = local_orbit_elements(position_km, velocity_km_s, mu_km3_s2)
%LOCAL_ORBIT_ELEMENTS Compute the stored Lambert orbit diagnostics.

position_km = position_km(:);
velocity_km_s = velocity_km_s(:);
specific_energy = 0.5 * dot(velocity_km_s, velocity_km_s) - ...
    mu_km3_s2 / norm(position_km);
semi_major_axis_km = -mu_km3_s2 / (2 * specific_energy);
eccentricity_vector = ((dot(velocity_km_s, velocity_km_s) - ...
    mu_km3_s2 / norm(position_km)) * position_km - ...
    dot(position_km, velocity_km_s) * velocity_km_s) / mu_km3_s2;
eccentricity = norm(eccentricity_vector);
end
