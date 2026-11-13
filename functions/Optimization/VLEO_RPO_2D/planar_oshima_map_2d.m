function [dv_vec, dDVdU] = planar_oshima_map_2d(U, is_case2b)
%PLANAR_OSHIMA_MAP_2D Map planar regularized impulse variables and their Jacobian.
%
% PROTOTYPE:
%   [dv_vec, dDVdU] = planar_oshima_map_2d(U, is_case2b)
%
% DESCRIPTION:
%   Use the nominal map [u^2-w^2; 2*u*w], or swap its components for the
%   case-2b branch. Either branch satisfies norm(delta-v) = u^2+w^2.
%
% INPUT:
%   U          [2 elements] Normalized regularized coefficients [u, w].
%   is_case2b  Logical selecting the alternate mapping branch.
%
% OUTPUT:
%   dv_vec     [2 x 1] Rotating-frame impulse [delta-vx; delta-vy] in VU.
%   dDVdU      [2 x 2] Analytic derivative of dv_vec with respect to [u, w].

U = U(:);
u = U(1);
w = U(2);

% Use the selected Oshima branch.
if is_case2b
    % Map the case-2b branch into planar delta-v.
    dv_vec = [2 * u * w; u^2 - w^2];

    dDVdU = [2 * w, 2 * u; 2 * u, -2 * w];
else
    % Map the nominal branch into planar delta-v.
    dv_vec = [u^2 - w^2; 2 * u * w];

    dDVdU = [2 * u, -2 * w; 2 * w, 2 * u];
end

end
