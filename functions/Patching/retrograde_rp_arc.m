function arc = retrograde_rp_arc(vinf_plus, rp, muE, rM_vec, vM_vec, sweep_rad, N)
%RETROGRADE_RP_ARC Sample a retrograde periapsis arc on the small circle.
%
% INPUTS
%   vinf_plus : 1x3 or 3x1 Earth-centered outgoing vInf+ vector
%   rp        : scalar Earth periapsis radius
%   muE       : scalar Earth GM
%   rM_vec    : 1x3 or 3x1 Moon position vector (sets rM=|rM_vec|)
%   vM_vec    : 1x3 or 3x1 Moon velocity vector (defines orbital plane)
%   sweep_rad : scalar total arc to sample on the small circle [rad]
%   N         : scalar number of samples
%
% OUTPUT
%   arc struct with fields:
%     Rp, Vp : Nx3 periapsis position/velocity (Earth-centered)
%

% Basic validation
if ~(isscalar(rp) && rp > 0)
    error('rp must be a positive scalar.');
end
if ~(isscalar(muE) && muE > 0)
    error('muE must be a positive scalar.');
end
if ~(isscalar(sweep_rad) && sweep_rad >= 0)
    error('sweep_rad must be a non-negative scalar.');
end
if ~(isscalar(N) && N > 0 && mod(N,1) == 0)
    error('N must be a positive integer.');
end

vinf_plus = vinf_plus(:).';
rM_vec = rM_vec(:).';
vM_vec = vM_vec(:).';

% Find Moon-circle intersections and pick retrograde center
sol = earthflyby_pick_zeta_intersect_moon_circle(vinf_plus, rp, muE, rM_vec, vM_vec);
if isempty(sol.zeta)
    arc = [];
    return
end

% Choose retrograde center (last solution)
zeta0 = sol.zeta(end);

% Build azimuth samples around zeta0
half = 0.5 * sweep_rad;
zetas = zeta0 + linspace(-half, +half, N);
zetas = wrap_to_pi(zetas);

% Sample periapsis states
Rp = zeros(N,3);
Vp = zeros(N,3);
for i = 1:N
    z = zetas(i);
    [Rp(i,:), Vp(i,:)] = rp_circle_from_vinfplus(vinf_plus, rp, muE, 1, z);
end

arc = struct('Rp', Rp, 'Vp', Vp);

end

function ang = wrap_to_pi(ang)
%WRAP_TO_PI Wrap azimuth offsets into the interval [-pi, pi).

    ang = mod(ang + pi, 2*pi) - pi;
end

