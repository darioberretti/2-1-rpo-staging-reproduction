function dxdt = rhs_SE3BP_3D(~, xx, mu_se)
%RHS_SE3BP_3D Evaluate normalized SE rotating state dynamics.
%
% PROTOTYPE:
%   dxdt = rhs_SE3BP_3D(t, xx, mu_se)
%   dxdt = rhs_SE3BP_3D(t, xx)
%
% DESCRIPTION:
%   Use autonomous CR3BP equations with primaries at x=-mu_se
%   and x=1-mu_se. Positions use DU_SE, velocity VU_SE=DU_SE/TU_SE
%   and derivatives use normalized time TU_SE. Run 04 supplies its
%   historical combined-mass parameter explicitly.
%
% INPUT:
%   t      Scalar integration time [TU_SE], unused.
%   xx     [6 elements] State [x;y;z;vx;vy;vz].
%   mu_se  Optional scalar secondary/total mass ratio [-]; omitted
%            or empty selects 1/(1+3.28900541e5).
%
% OUTPUT:
%   dxdt  [6 x 1] Derivative in input state order.
%            Equations/Jacobian are evaluated by rhs_SE3BP_3D_fA.

    if nargin < 3 || isempty(mu_se)
        mu_se = default_mu_se();
    end

    validateattributes(xx, {'numeric'}, {'real', 'finite', 'vector', 'numel', 6}, ...
        mfilename, 'xx', 2);
    validateattributes(mu_se, {'numeric'}, {'real', 'finite', 'scalar', '>', 0, '<', 0.5}, ...
        mfilename, 'mu_se', 3);

    [dxdt, ~] = rhs_SE3BP_3D_fA(xx(:), mu_se);
end

function mu_se = default_mu_se()
%DEFAULT_MU_SE Return the mass ratio implied by the shared solar mass scaling.

    MU_SCALED_SUN = 3.28900541e5;
    mu_se = 1 / (1 + MU_SCALED_SUN);
end
