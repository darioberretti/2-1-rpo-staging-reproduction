function interplanetaryEarthMoon = interplanetary_MoonEarth_patch_opt(vInfPlusEarth, N, hMinEarthFlyby, sweep_rad, type, nThetaMoon)
%INTERPLANETARY_MOONEARTH_PATCH_OPT Build Moon-Earth patched-conic candidates.
%
% PROTOTYPE:
%   database = interplanetary_MoonEarth_patch_opt(vInfPlusEarth, N, ...
%       hMinEarthFlyby, sweep_rad, type, nThetaMoon)
%
% DESCRIPTION:
%   Sample Earth-perigee directions and Moon orbital phases, then find a
%   Lambert flight time whose terminal velocity is perpendicular to radius.
%   Build lunar incoming/outgoing matching quantities for the patch selector.
%
% INPUT:
%   vInfPlusEarth    [3 elements] Outgoing Earth-relative ECI v-infinity [km/s].
%   N                Number of Earth-perigee samples.
%   hMinEarthFlyby   Earth-perigee altitude [km].
%   sweep_rad        Perigee arc span [rad] for retrograde sampling.
%   type             0 for a full circle, 1 for a retrograde arc.
%   nThetaMoon       Optional Moon-phase sample count, default 100.
%
% OUTPUT:
%   interplanetaryEarthMoon  Table of valid patch candidates. R1/R2/Rp use km,
%                             VI/VF/Vp and v-infinity vectors use km/s, ToF
%                             uses seconds, and ToF_days uses days. Angles
%                             retain the explicitly named source conventions.

if nargin < 6 || isempty(nThetaMoon)
    nThetaMoon = 100;
end
% Constants and setup
rEM = 384400;               % [km] Earth-Moon distance
muE = 3.98600433e+05;       % [km^3/s^2] Earth GM
RE  = astroConstants(23);   % [km] Earth radius
h   = hMinEarthFlyby;       % [km]
rp  = RE + h;               % [km]
day = 3600*24;              % [s]

vinf_plus = vInfPlusEarth(:).';
% Periapsis sampling
if ~type
    [Rp, Vp] = rp_circle_from_vinfplus(vinf_plus, rp, muE, N);
else
    rM_vec = [rEM 0 0]';
    vM_vec = [0 sqrt(muE/rEM) 0]';
    arc = retrograde_rp_arc(vinf_plus, rp, muE, rM_vec, vM_vec, sweep_rad, N);
    if ~isempty(arc)
        Rp = arc.Rp;
        Vp = arc.Vp;
    else
        [Rp, Vp] = rp_circle_from_vinfplus(vinf_plus, rp, muE, N);
    end
end

% Initial ToF scale (heuristic)
TOF0 = pi * sqrt((rp + rEM)^3 / muE);  % [s]
% Search settings
% tol is normalized dot-residual on perpendicularity:
% err = |dot(VF,R2)| / (||VF||*||R2||)

tol = 1e-3;
maxIterTOF = 60;

% Moon inertial velocity (constant in patched conics frame)
vMoon = [0 sqrt(muE/rEM) 0];
% Preallocate outputs
M = N * nThetaMoon;
R1_all          = zeros(M,3);
R2_all          = zeros(M,3);
VI_all          = zeros(M,3);
VF_all          = zeros(M,3);
ToF_all         = zeros(M,1);
ToF_days_all    = zeros(M,1);
vInfMoon_all    = zeros(M,3);
vInfMoonMag_all = zeros(M,1);
a_all           = zeros(M,1);
e_all           = zeros(M,1);
eps_all         = zeros(M,1);
theta_all       = zeros(M,1);
dvp_all         = zeros(M,3);
dvpMag_all      = zeros(M,1);
Vp_all          = zeros(M,3);
% Loop over periapsis samples
for j = 1:N
    R2 = Rp(j,:);
    vp = Vp(j,:);
    r2 = norm(R2);

    % Moon phase sweep (aligned with R2 azimuth)
    azimuth_R2 = atan2(R2(2), R2(1));
    thetaMoon0 = wrapTo2Pi(azimuth_R2 - pi);
    if thetaMoon0 > 3/2*pi && thetaMoon0 <= 2*pi
        thetaMoon0 = thetaMoon0 - 2*pi;
    end
    thetaMoonf = thetaMoon0 - pi/2;

    thetaMoonVect = linspace(thetaMoon0, thetaMoonf, nThetaMoon+1);
    thetaMoonVect(end) = [];

    cth = cos(thetaMoonVect);
    sth = sin(thetaMoonVect);

    % Precompute all R1 and chord-based ToF seeds
    R1s = rEM * [cth(:), sth(:), zeros(nThetaMoon,1)];
    chordT = vecnorm(R1s - R2, 2, 2) ./ sqrt(muE);
    chordT(chordT < 1) = 1;

    % Global write range for this block
    base = (j-1)*nThetaMoon;
    idxRange = base + (1:nThetaMoon);

    % Per-block arrays
    R1_blk          = zeros(nThetaMoon,3);
    R2_blk          = zeros(nThetaMoon,3);
    VI_blk          = zeros(nThetaMoon,3);
    VF_blk          = zeros(nThetaMoon,3);
    ToF_blk         = zeros(nThetaMoon,1);
    ToF_days_blk    = zeros(nThetaMoon,1);
    vInfMoon_blk    = zeros(nThetaMoon,3);
    vInfMoonMag_blk = zeros(nThetaMoon,1);
    a_blk           = zeros(nThetaMoon,1);
    e_blk           = zeros(nThetaMoon,1);
    eps_blk         = zeros(nThetaMoon,1);
    theta_blk       = zeros(nThetaMoon,1);
    dvp_blk         = zeros(nThetaMoon,3);
    dvpMag_blk      = zeros(nThetaMoon,1);
    Vp_blk          = zeros(nThetaMoon,3);

    for k = 1:nThetaMoon
        thetaMoon = thetaMoonVect(k);
        R1 = R1s(k,:);

        % Bracketed root-find on g(ToF) = dot(VF,R2) = 0
        [A, E, VI, VF, ToF_star] = RPO2DS_solve_perp_lambert_tof(R1, R2, muE, chordT(k), TOF0, r2, tol, maxIterTOF);

        % Specific orbital energy at |R1|=rEM
        eps = 0.5*dot(VI,VI) - muE/rEM;

        % vInfMoon via local R/T/H at R1
        rdir = R1 / norm(R1);
        hdir = [0 0 1]';
        tdir = cross(hdir, rdir);
        VI_RTH = [dot(VI',rdir), dot(VI',tdir), dot(VI',hdir)];
        vInfMoon    = VI_RTH(:).' - vMoon;
        vInfMoonMag = norm(vInfMoon);

        % Maneuver at pericenter
        deltavp    = - VF(:).' + vp(:).';
        deltavpMag = norm(deltavp);

        % Store
        R1_blk(k,:)          = R1;
        R2_blk(k,:)          = R2;
        VI_blk(k,:)          = VI;
        VF_blk(k,:)          = VF;
        ToF_blk(k)           = ToF_star;
        ToF_days_blk(k)      = ToF_star / day;
        vInfMoon_blk(k,:)    = vInfMoon;
        vInfMoonMag_blk(k)   = vInfMoonMag;
        a_blk(k)             = A;
        e_blk(k)             = E;
        eps_blk(k)           = eps;
        theta_blk(k)         = thetaMoon;
        dvp_blk(k,:)         = deltavp;
        dvpMag_blk(k)        = deltavpMag;
        Vp_blk(k,:)          = vp;
    end

    % Copy block results into global arrays
    R1_all(idxRange,:)          = R1_blk;
    R2_all(idxRange,:)          = R2_blk;
    VI_all(idxRange,:)          = VI_blk;
    VF_all(idxRange,:)          = VF_blk;
    ToF_all(idxRange)           = ToF_blk;
    ToF_days_all(idxRange)      = ToF_days_blk;
    vInfMoon_all(idxRange,:)    = vInfMoon_blk;
    vInfMoonMag_all(idxRange)   = vInfMoonMag_blk;
    a_all(idxRange)             = a_blk;
    e_all(idxRange)             = e_blk;
    eps_all(idxRange)           = eps_blk;
    theta_all(idxRange)         = theta_blk;
    dvp_all(idxRange,:)         = dvp_blk;
    dvpMag_all(idxRange)        = dvpMag_blk;
    Vp_all(idxRange,:)          = Vp_blk;
end
% Build output table
MoonEarthDatabase = table( ...
    R1_all, R2_all, VI_all, VF_all, ToF_all, ToF_days_all, ...
    vInfMoon_all, vInfMoonMag_all, ...
    a_all, e_all, eps_all, theta_all, dvp_all, dvpMag_all, Vp_all, ...
    'VariableNames', {'R1','R2','VI','VF','ToF','ToF_days', ...
                      'vInfMoon','vInfMoonMag','a','e','eps','thetaMoon', ...
                      'deltavp','deltavpMag','Vp'} );

% Zero-out tiny values in selected columns
zeroTol = 1e-10;
for v = intersect(["R1","R2","VI","VF","ToF","vInfMoon"], string(MoonEarthDatabase.Properties.VariableNames))
    MoonEarthDatabase.(v)(abs(MoonEarthDatabase.(v)) < zeroTol) = 0;
end

interplanetaryEarthMoon = MoonEarthDatabase;

end
