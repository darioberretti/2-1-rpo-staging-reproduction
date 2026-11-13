function [xf, STMf, dxf_dtheta] = propagate_bcr4bp_2d_stm_theta_fast(x0, t0, tf, mu, thetaSun0, absTol, relTol, compute_jacobian)
%PROPAGATE_BCR4BP_2D_STM_THETA_FAST Propagate a planar state with optional STM and Sun-phase sensitivity.
%
% PROTOTYPE:
%   [xf, STMf, dxf_dtheta] = propagate_bcr4bp_2d_stm_theta_fast( ...
%       x0, t0, tf, mu, thetaSun0, absTol, relTol, compute_jacobian)
%
% DESCRIPTION:
%   Integrate the BCR4BP state, identity-initialized STM, and initially zero
%   derivative with respect to thetaSun0 in one augmented pass. The STM maps
%   segment-start state perturbations to its end. Negligible time intervals
%   return the supplied state and initial sensitivity values directly.
%
% INPUT:
%   x0                [4 elements] Rotating barycentric [x, y, vx, vy] in DU/VU.
%   t0, tf            Segment-start/final epochs [TU].
%   mu                Secondary-to-total Earth-Moon mass ratio [-].
%   thetaSun0         Initial solar angle [rad]; thetaSun = thetaSun0 + ws*t.
%   absTol, relTol    ode113 integration tolerances.
%   compute_jacobian  Optional logical, default true; false integrates only state.
%
% OUTPUT:
%   xf                [4 x 1] Final state in DU/VU.
%   STMf              [4 x 4] Final state-transition matrix, or [] when disabled.
%   dxf_dtheta        [4 x 1] Final-state derivative per radian of thetaSun0,
%                       or [] when disabled.

if nargin < 8 || isempty(compute_jacobian)
    compute_jacobian = true;
end

% Reshape the initial state into a column vector.
x0 = x0(:);

% Return the initial values directly when the segment duration is zero.
if abs(tf - t0) < 1e-12
    xf = x0;
    if compute_jacobian
        STMf = eye(4);
        dxf_dtheta = zeros(4, 1);
    else
        STMf = [];
        dxf_dtheta = [];
    end
    return;
end

% Build the ODE options once for the selected propagation mode.
ode_opts = odeset('RelTol', relTol, 'AbsTol', absTol);

% Integrate the augmented system only when the Jacobian is needed.
if compute_jacobian
    % Append the identity STM and the zero theta sensitivity to the state.
    Phi0 = eye(4);
    z0 = [x0; Phi0(:); zeros(4, 1)];

    % Integrate the full augmented system over the current segment.
    [~, zz] = ode113(@(t, z) rhs_BCR4BP_2D_STM_theta_fast(t, z, mu, thetaSun0), [t0, tf], z0, ode_opts);

    % Read the propagated state, STM, and theta sensitivity.
    xf = zz(end, 1:4).';
    STMf = reshape(zz(end, 5:20), 4, 4);
    dxf_dtheta = zz(end, 21:24).';
else
    % Integrate only the physical state when the Jacobian is not needed.
    [~, xx] = ode113(@(t, x) rhs_BCR4BP_2D(t, x, mu, thetaSun0), [t0, tf], x0, ode_opts);

    % Read the propagated state and leave the sensitivity outputs empty.
    xf = xx(end, 1:4).';
    STMf = [];
    dxf_dtheta = [];
end

end
