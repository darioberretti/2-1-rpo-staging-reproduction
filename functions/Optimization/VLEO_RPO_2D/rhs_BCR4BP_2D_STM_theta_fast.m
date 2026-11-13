function dzdt = rhs_BCR4BP_2D_STM_theta_fast(t, z, mu, thetaSun0)
%RHS_BCR4BP_2D_STM_THETA_FAST Evaluate planar state, STM, and initial-Sun-angle sensitivity dynamics.
%
% PROTOTYPE:
%   dzdt = rhs_BCR4BP_2D_STM_theta_fast(t, z, mu, thetaSun0)
%
% DESCRIPTION:
%   Use the frozen Earth-Moon rotating BCR4BP. Alongside dx/dt = f,
%   integrate dPhi/dt = A*Phi and ds/dt = A*s + df/dthetaSun0, with
%   A = df/dx. The flattened STM uses MATLAB column-major order. These
%   are variational equations; no extra full-mesh finite difference is used.
%
% INPUT:
%   t          Current normalized epoch [TU].
%   z          [24 x 1] [state(4); Phi(:)(16); s(4)], where state is
%                [x; y; vx; vy] in DU/VU and s = dstate/dthetaSun0.
%   mu         Earth-Moon mass ratio [-].
%   thetaSun0  Initial Sun phase [rad]; solar angle is thetaSun0 + ws*t.
%
% OUTPUT:
%   dzdt       [24 x 1] Augmented derivative in the same order, per TU.

x = z(1);
y = z(2);
vx = z(3);
vy = z(4);
Phi = reshape(z(5:20), 4, 4);

% Read the propagated sensitivity with respect to thetaSun0.
dx_dtheta = z(21:24);

% Store the fixed Sun parameters used by the normalized BCR4BP model.
ms = 3.28900541e5;
rho = 3.88811143e2;
ws = -9.25195985e-1;

% Build the absolute Sun phase at the current epoch.
thetaSun = thetaSun0 + ws * t;
cSun = cos(thetaSun);
sSun = sin(thetaSun);

% Build the body-relative distances.
r1_sq = (x + mu)^2 + y^2;
r2_sq = (x + mu - 1)^2 + y^2;
dx3 = x - rho * cSun;
dy3 = y - rho * sSun;
r3_sq = dx3^2 + dy3^2;

% Build the inverse powers used by the accelerations and Jacobians.
r1 = sqrt(r1_sq);
r2 = sqrt(r2_sq);
r3 = sqrt(r3_sq);
r1_inv3 = 1 / (r1_sq * r1);
r2_inv3 = 1 / (r2_sq * r2);
r3_inv3 = 1 / (r3_sq * r3);
r1_inv5 = r1_inv3 / r1_sq;
r2_inv5 = r2_inv3 / r2_sq;
r3_inv5 = r3_inv3 / r3_sq;

% Build the state derivatives of the effective potential.
dUdx = x ...
    - (1 - mu) * (mu + x) * r1_inv3 ...
    + mu * (1 - mu - x) * r2_inv3 ...
    - ms * dx3 * r3_inv3 ...
    - ms / rho^2 * cSun;
dUdy = y ...
    - (1 - mu) * y * r1_inv3 ...
    - mu * y * r2_inv3 ...
    - ms * dy3 * r3_inv3 ...
    - ms / rho^2 * sSun;

% Assemble the 2D state dynamics.
state_dot = zeros(4, 1);
state_dot(1) = vx;
state_dot(2) = vy;
state_dot(3) = dUdx + 2 * vy;
state_dot(4) = dUdy - 2 * vx;

% Build the 2D state Jacobian A = df/dx.
A = zeros(4, 4);
A(1, 3) = 1;
A(2, 4) = 1;
A(3, 1) = 1 ...
    - (1 - mu) * r1_inv3 + 3 * (1 - mu) * (x + mu)^2 * r1_inv5 ...
    - mu * r2_inv3 + 3 * mu * (x + mu - 1)^2 * r2_inv5 ...
    - ms * r3_inv3 + 3 * ms * dx3^2 * r3_inv5;
A(3, 2) = 3 * (1 - mu) * (x + mu) * y * r1_inv5 ...
    + 3 * mu * (x + mu - 1) * y * r2_inv5 ...
    + 3 * ms * dx3 * dy3 * r3_inv5;
A(3, 4) = 2;
A(4, 1) = A(3, 2);
A(4, 2) = 1 ...
    - (1 - mu) * r1_inv3 + 3 * (1 - mu) * y^2 * r1_inv5 ...
    - mu * r2_inv3 + 3 * mu * y^2 * r2_inv5 ...
    - ms * r3_inv3 + 3 * ms * dy3^2 * r3_inv5;
A(4, 3) = -2;

% Build the direct thetaSun0 sensitivity df/dthetaSun0.
dr3sq_dtheta = 2 * rho * (dx3 * sSun - dy3 * cSun);
dr3inv3_dtheta = -1.5 * r3_inv5 * dr3sq_dtheta;
ddx3_dtheta = rho * sSun;
ddy3_dtheta = -rho * cSun;

% Differentiate the Sun contribution in the x acceleration.
ddUdx_dtheta = -ms * (ddx3_dtheta * r3_inv3 + dx3 * dr3inv3_dtheta) ...
    + ms / rho^2 * sSun;

% Differentiate the Sun contribution in the y acceleration.
ddUdy_dtheta = -ms * (ddy3_dtheta * r3_inv3 + dy3 * dr3inv3_dtheta) ...
    - ms / rho^2 * cSun;

% Assemble the direct thetaSun0 forcing vector.
df_dtheta = zeros(4, 1);
df_dtheta(3) = ddUdx_dtheta;
df_dtheta(4) = ddUdy_dtheta;

% Differentiate the STM and the propagated theta sensitivity.
Phi_dot = A * Phi;
dx_dtheta_dot = A * dx_dtheta + df_dtheta;

% Assemble the full augmented derivative.
dzdt = zeros(24, 1);
dzdt(1:4) = state_dot;
dzdt(5:20) = Phi_dot(:);
dzdt(21:24) = dx_dtheta_dot;

end
