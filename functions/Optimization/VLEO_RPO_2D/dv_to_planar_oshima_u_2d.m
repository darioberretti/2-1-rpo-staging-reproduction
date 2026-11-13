function [U, is_case2b] = dv_to_planar_oshima_u_2d(dv_vec)
%DV_TO_PLANAR_OSHIMA_U_2D Invert one planar impulse into an Oshima warm-start branch.
%
% PROTOTYPE:
%   [U, is_case2b] = dv_to_planar_oshima_u_2d(dv_vec)
%
% DESCRIPTION:
%   Choose the nominal inverse except near the pure negative-x impulse,
%   where case 2b avoids the singular nominal denominator. For negligible
%   impulses retain the implemented nonzero eps seed.
%
% INPUT:
%   dv_vec     [2 elements] Rotating-frame impulse [delta-vx, delta-vy] in VU.
%
% OUTPUT:
%   U          [2 x 1] Regularized coefficients for planar_oshima_map_2d.
%   is_case2b  Logical branch flag; keep it fixed with this warm start.

dv_vec = dv_vec(:);
dvx = dv_vec(1);
dvy = dv_vec(2);

% Build the useful planar magnitudes.
dv_mag = hypot(dvx, dvy);
tol = 1e-14;
is_case2b = false;

% Use the zero-delta-v branch when needed.
if dv_mag < tol
    % Keep a numerically nonzero regularized seed.
    u = eps;
    w = eps;

% Use the nominal branch when dvx is not the pure negative-axis case.
elseif dvx > 0 || abs(dvy) > tol
    % Recover the nominal branch from the planar magnitude.
    u = sqrt(max(eps, 0.5 * (dv_mag + dvx)));
    w = dvy / (2 * u);

% Use the case-2b branch for a pure negative-x impulse.
else
    % Recover the case-2b branch directly.
    u = sqrt(max(eps, -0.5 * dvx));
    w = -sqrt(max(eps, -0.5 * dvx));
    is_case2b = true;
end

U = [u; w];

end
