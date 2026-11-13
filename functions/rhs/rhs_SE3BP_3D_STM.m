function dxdt = rhs_SE3BP_3D_STM(~, xx, mu_se)
%RHS_SE3BP_3D_STM Evaluate normalized SE rotating state and STM dynamics.
%
% PROTOTYPE:
%   dxdt = rhs_SE3BP_3D_STM(t, xx, mu_se)
%   dxdt = rhs_SE3BP_3D_STM(t, xx)
%
% DESCRIPTION:
%   Use autonomous CR3BP equations with primaries at x=-mu_se
%   and x=1-mu_se. Positions use DU_SE, velocity VU_SE=DU_SE/TU_SE
%   and derivatives use normalized time TU_SE. Run 04 supplies its
%   historical combined-mass parameter explicitly.
%   Phi_dot=A*Phi; an identity-initialized Phi maps initial-state
%   perturbations to current-state perturbations. Flatten Phi in
%   MATLAB column-major order.
%
% INPUT:
%   t      Scalar integration time [TU_SE], unused.
%   xx     [42 elements] Augmented [x;y;z;vx;vy;vz;Phi(:)], Phi [6 x 6].
%   mu_se  Optional scalar secondary/total mass ratio [-]; omitted
%            or empty selects 1/(1+3.28900541e5).
%
% OUTPUT:
%   dxdt  [42 x 1] Derivative [state_dot;Phi_dot(:)].
%            Equations/Jacobian are evaluated by rhs_SE3BP_3D_fA.

    if nargin < 3 || isempty(mu_se)
        mu_se = default_mu_se();
    end

    validateattributes(xx, {'numeric'}, {'real', 'finite', 'vector', 'numel', 42}, ...
        mfilename, 'xx', 2);
    validateattributes(mu_se, {'numeric'}, {'real', 'finite', 'scalar', '>', 0, '<', 0.5}, ...
        mfilename, 'mu_se', 3);

    % Unpack the state/STM and apply the variational equation.
    xx = xx(:);
    x_state = xx(1:6);
    Phi = reshape(xx(7:end), 6, 6);

    [f_state, A] = rhs_SE3BP_3D_fA(x_state, mu_se);
    Phi_dot = A * Phi;

    dxdt = [f_state; Phi_dot(:)];
end

function mu_se = default_mu_se()
%DEFAULT_MU_SE Return the mass ratio implied by the shared solar mass scaling.

    MU_SCALED_SUN = 3.28900541e5;
    mu_se = 1 / (1 + MU_SCALED_SUN);
end
