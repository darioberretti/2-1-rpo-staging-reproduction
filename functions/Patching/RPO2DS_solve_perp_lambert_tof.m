function [A, E, VI, VF, ToFStar] = RPO2DS_solve_perp_lambert_tof(R1, R2, mu, ToFSeed, TOFscale, r2norm, tolLocal, maxIt, maxTofSec)
%RPO2DS_SOLVE_PERP_LAMBERT_TOF Find a Lambert time with terminal velocity normal to radius.
%
% PROTOTYPE:
%   [A, E, VI, VF, ToFStar] = RPO2DS_solve_perp_lambert_tof(R1, R2, mu, ...
%       ToFSeed, TOFscale, r2norm, tolLocal, maxIt, maxTofSec)
%
% DESCRIPTION:
%   Bracket a zero of dot(VF, R2) by expanding toward longer/shorter times,
%   then refine it while preserving the chosen Lambert orientation. Normalize
%   the perpendicularity residual by the radius and speed magnitudes.
%
% INPUT:
%   R1, R2     [3 elements] Endpoint positions [km], in the same inertial frame.
%   mu         Central-body gravitational parameter [km^3/s^2].
%   ToFSeed    Initial flight time [s].
%   TOFscale   Time scale used to expand the bracket [s].
%   r2norm     Optional norm(R2) [km].
%   tolLocal   Optional normalized angular residual tolerance, default 1e-3.
%   maxIt      Optional iteration limit, default 60.
%   maxTofSec  Optional maximum flight time [s], default Inf.
%
% OUTPUT:
%   A, E       Selected semi-major axis [km] and eccentricity [-].
%   VI, VF     Departure/arrival velocities [km/s].
%   ToFStar    Selected Lambert flight time [s].

R1 = R1(:);
R2 = R2(:);
if nargin < 6 || isempty(r2norm)
    r2norm = norm(R2);
end

% Use standard local tolerances unless the caller requests tighter settings.
if nargin < 7 || isempty(tolLocal)
    tolLocal = 1e-3;
end
if nargin < 8 || isempty(maxIt)
    maxIt = 60;
end
if nargin < 9 || isempty(maxTofSec)
    maxTofSec = inf;
end

% Evaluate the initial Lambert seed.
Tlo = min(max(1, ToFSeed), maxTofSec);
[Alo, Elo, VIlo, VFlo] = local_lambert(R1, R2, Tlo, mu);
glo = dot(VFlo, R2);
errLo = abs(glo) / (norm(VFlo) * r2norm);
if isfinite(errLo) && errLo < tolLocal
    A = Alo; E = Elo; VI = VIlo; VF = VFlo; ToFStar = Tlo;
    return
end

% Expand the bracket toward longer transfer times.
step = 0.5 * TOFscale;
Thi = min(Tlo + step, maxTofSec);
[Ahi, Ehi, VIhi, VFhi] = local_lambert(R1, R2, Thi, mu);
ghi = dot(VFhi, R2);
flips = false;
for k = 1:8
    if sign(ghi) ~= sign(glo) && isfinite(ghi) && isfinite(glo)
        flips = true;
        break
    end
    step = step * 2;
    Thi = min(Tlo + step, maxTofSec);
    [Ahi, Ehi, VIhi, VFhi] = local_lambert(R1, R2, Thi, mu);
    ghi = dot(VFhi, R2);
    if Thi >= maxTofSec
        break
    end
end

% Try shorter transfer times when the longer-time bracket did not cross.
if ~flips
    step = 0.5 * TOFscale;
    Thi = max(1, Tlo - step);
    [Ahi, Ehi, VIhi, VFhi] = local_lambert(R1, R2, Thi, mu);
    ghi = dot(VFhi, R2);
    for k = 1:8
        if sign(ghi) ~= sign(glo) && isfinite(ghi) && isfinite(glo)
            flips = true;
            break
        end
        step = step * 2;
        Thi = max(1, Tlo - step);
        [Ahi, Ehi, VIhi, VFhi] = local_lambert(R1, R2, Thi, mu);
        ghi = dot(VFhi, R2);
    end
end

% Keep the lowest residual endpoint when no sign-changing bracket exists.
if ~flips
    errHi = abs(ghi) / (norm(VFhi) * r2norm);
    if errHi < errLo
        A = Ahi; E = Ehi; VI = VIhi; VF = VFhi; ToFStar = Thi;
    else
        A = Alo; E = Elo; VI = VIlo; VF = VFlo; ToFStar = Tlo;
    end
    return
end

% Bisect the bracket and retain the best endpoint state for the final row.
AloB = Alo; EloB = Elo; VIloB = VIlo; VFloB = VFlo;
AhiB = Ahi; EhiB = Ehi; VIhiB = VIhi; VFhiB = VFhi;
for it = 1:maxIt
    Tmid = 0.5 * (Tlo + Thi);
    [Amid, Emid, VImid, VFmid] = local_lambert(R1, R2, Tmid, mu);
    gmid = dot(VFmid, R2);
    err = abs(gmid) / (norm(VFmid) * r2norm);
    if ~isfinite(err)
        Thi = Tmid; AhiB = Amid; EhiB = Emid; VIhiB = VImid; VFhiB = VFmid;
        continue
    end
    if isfinite(err) && err < tolLocal
        A = Amid; E = Emid; VI = VImid; VF = VFmid; ToFStar = Tmid;
        return
    end
    if sign(gmid) == sign(glo)
        Tlo = Tmid; glo = gmid; AloB = Amid; EloB = Emid; VIloB = VImid; VFloB = VFmid;
    else
        Thi = Tmid; AhiB = Amid; EhiB = Emid; VIhiB = VImid; VFhiB = VFmid;
    end
end

% Use the best bracket endpoint after the iteration cap.
errLo = abs(dot(VFloB, R2)) / (norm(VFloB) * r2norm);
errHi = abs(dot(VFhiB, R2)) / (norm(VFhiB) * r2norm);
if errHi < errLo
    A = AhiB; E = EhiB; VI = VIhiB; VF = VFhiB; ToFStar = Thi;
else
    A = AloB; E = EloB; VI = VIloB; VF = VFloB; ToFStar = Tlo;
end
end

function [semi_major_axis_km, eccentricity, departure_velocity_km_s, arrival_velocity_km_s] = ...
        local_lambert(R1, R2, time_of_flight_seconds, mu_km3_s2)
%LOCAL_LAMBERT Use the redistributable solver for the historical retrograde leg.

R1 = R1(:);
R2 = R2(:);

[departure_velocity_km_s, arrival_velocity_km_s] = rpo2ds.solve_lambert( ...
    R1, R2, time_of_flight_seconds / 86400, mu_km3_s2, "retrograde");
radius_km = norm(R1);
specific_energy = 0.5 * dot(departure_velocity_km_s, departure_velocity_km_s) ...
    - mu_km3_s2 / radius_km;
semi_major_axis_km = -mu_km3_s2 / (2 * specific_energy);
eccentricity_vector = ((dot(departure_velocity_km_s, departure_velocity_km_s) ...
    - mu_km3_s2 / radius_km) * R1 - ...
    dot(R1, departure_velocity_km_s) * departure_velocity_km_s) / mu_km3_s2;
eccentricity = norm(eccentricity_vector);
end
