function [dr1dot_dt, dr1dot_dx, dr1dot_dy, dr1dot_dz, ...
          dr1dot_dvx, dr1dot_dvy, dr1dot_dvz, ...
          dr1ddot_dt, dr1ddot_dx, dr1ddot_dy, dr1ddot_dz, ...
          dr1ddot_dvx, dr1ddot_dvy, dr1ddot_dvz] = get_dr1ddot_dtX(X, constants)
%GET_DR1DDOT_DTX Differentiate primary-1 radial velocity and acceleration.
%
% PROTOTYPE:
%   [dr1dot_dt, dr1dot_dx, dr1dot_dy, dr1dot_dz, ...
%       dr1dot_dvx, dr1dot_dvy, dr1dot_dvz, ...
%       dr1ddot_dt, dr1ddot_dx, dr1ddot_dy, dr1ddot_dz, ...
%       dr1ddot_dvx, dr1ddot_dvy, dr1ddot_dvz] = get_dr1ddot_dtX(X, constants)
%
% DESCRIPTION:
%   Evaluate analytic partials of dr1/dt and d^2r1/dt^2 with respect
%   to explicit time and the six normalized rotating-frame state components.
%   Explicit time derivatives are zero for the autonomous CR3BP.
%
% INPUT:
%   X          [6 elements] State [x, y, z, vx, vy, vz] in DU/VU.
%   constants  Legacy vector whose first entry is the mass ratio mu [-].
%
% OUTPUT:
%   First 7 scalars   Partials of r1dot with respect to
%                       [t, x, y, z, vx, vy, vz], in that order.
%   Last 7 scalars    Corresponding partials of r1ddot. Each derivative
%                       uses the output/input units from get_r1ddot.

mu = constants(1);

x  = X(1);
y  = X(2);
z  = X(3);

vx = X(4);
vy = X(5);
vz = X(6);

% Compute acceleration from CR3BP equations of motion
% dX = [vx; vy; vz; ax; ay; az]
dX = rhs_CR3BP_3D(0, X, mu);
ax = dX(4);
ay = dX(5);
az = dX(6);

% Useful radius terms (distances to primaries)
% Primary 1 at x = -mu
% Primary 2 at x = 1 - mu
dx1 = x + mu;        % displacement w.r.t. primary 1
dx2 = x - 1 + mu;    % displacement w.r.t. primary 2

% Distance factors for primary 1 (the one we care about for r1, r1dot, r1ddot)
r1_05 = (dx1^2 + y^2 + z^2)^(-0.5);   % 1 / |r1|
r1_15 = (dx1^2 + y^2 + z^2)^(-1.5);   % 1 / |r1|^3
r1_25 = (dx1^2 + y^2 + z^2)^(-2.5);   % 1 / |r1|^5

% Distance factors for primary 2 (still needed for second derivatives of U)
r2_15 = (dx2^2 + y^2 + z^2)^(-1.5);   % 1 / |r2|^3
r2_25 = (dx2^2 + y^2 + z^2)^(-2.5);   % 1 / |r2|^5

% Second derivatives of the effective potential U
% Uxx, Uxy, Uxz, ...
Uxx = -1 ...
      - 3*mu*(dx2)^2 * r2_25 + mu * r2_15 ...
      - 3*(1 - mu)*(dx1)^2 * r1_25 + (1 - mu)*r1_15;

Uxy = -3*mu*dx2*y * r2_25 ...
      -3*(1 - mu)*dx1*y * r1_25;

Uxz = -3*mu*dx2*z * r2_25 ...
      -3*(1 - mu)*dx1*z * r1_25;

Uyx = Uxy;

Uyy = -1 ...
      -3*mu*y^2 * r2_25 + mu * r2_15 ...
      -3*(1 - mu)*y^2 * r1_25 + (1 - mu)*r1_15;

Uyz = -3*(1 - mu)*y*z * r1_25 ...
      -3*mu*y*z * r2_25;

Uzx = Uxz;
Uzy = Uyz;

Uzz = -3*(1 - mu)*z^2 * r1_25 ...
      -3*mu*z^2 * r2_25 ...
      + (1 - mu)*r1_15 + mu*r2_15;

% First derivatives of r1dot

% r1dot does not explicitly depend on time
dr1dot_dt  = 0;

% ∂(r1dot)/∂x, ∂(r1dot)/∂y, ∂(r1dot)/∂z
dr1dot_dx = r1_15 * ( vx*(y^2 + z^2)   - dx1*(y*vy + z*vz) );
dr1dot_dy = r1_15 * ( vy*(dx1^2 + z^2) - y*(dx1*vx + z*vz) );
dr1dot_dz = r1_15 * ( vz*(dx1^2 + y^2) - z*(dx1*vx + y*vy) );

% ∂(r1dot)/∂vx, ∂(r1dot)/∂vy, ∂(r1dot)/∂vz
dr1dot_dvx = dx1 * r1_05;
dr1dot_dvy = y   * r1_05;
dr1dot_dvz = z   * r1_05;

% Second derivatives needed for r1ddot derivatives

% r1dot_xx, r1dot_xy, r1dot_xz, ...
r2dot_xx = -3 * r1_25 * dx1 * ( vx*(y^2 + z^2) - dx1*(y*vy + z*vz) ) ...
           - r1_15 * (y*vy + z*vz);

r2dot_xy = -3 * r1_25 * y * ( vx*(y^2 + z^2) - dx1*(y*vy + z*vz) ) ...
           + r1_15 * (2*y*vx - dx1*vy);

r2dot_xz = -3 * r1_25 * z * ( vx*(y^2 + z^2) - dx1*(y*vy + z*vz) ) ...
           + r1_15 * (2*z*vx - dx1*vz);

r2dot_yx = r2dot_xy;

r2dot_yy = -3 * r1_25 * y * ( vy*(dx1^2 + z^2) - y*(dx1*vx + z*vz) ) ...
           - r1_15 * (dx1*vx + z*vz);

r2dot_yz = -3 * r1_25 * z * ( vy*(dx1^2 + z^2) - y*(dx1*vx + z*vz) ) ...
           + r1_15 * (2*z*vy - y*vz);

r2dot_zx = r2dot_xz;
r2dot_zy = r2dot_yz;

r2dot_zz = -3 * r1_25 * z * ( vz*(dx1^2 + y^2) - z*(dx1*vx + y*vy) ) ...
           - r1_15 * (dx1*vx + y*vy);

% Mixed derivatives with velocities
r2dot_xvx =  r1_15 * (y^2 + z^2);
r2dot_xvy = -r1_15 * dx1 * y;
r2dot_xvz = -r1_15 * dx1 * z;

r2dot_yvx = -r1_15 * dx1 * y;
r2dot_yvy =  r1_15 * (dx1^2 + z^2);
r2dot_yvz = -r1_15 * y * z;

r2dot_zvx = -r1_15 * dx1 * z;
r2dot_zvy = -r1_15 * y * z;
r2dot_zvz =  r1_15 * (dx1^2 + y^2);

% Derivatives of r1ddot

% r1ddot does not explicitly depend on time in this formulation
dr1ddot_dt = 0;

% ∂(r1ddot)/∂x
dr1ddot_dx = r2dot_xx * vx + r2dot_xy * vy + r2dot_xz * vz ...
             + r2dot_xvx * ax - dr1dot_dvx * Uxx ...
             + r2dot_xvy * ay - dr1dot_dvy * Uyx ...
             + r2dot_xvz * az - dr1dot_dvz * Uzx;

% ∂(r1ddot)/∂y
dr1ddot_dy = r2dot_yx * vx + r2dot_yy * vy + r2dot_yz * vz ...
             + r2dot_yvx * ax - dr1dot_dvx * Uyx ...
             + r2dot_yvy * ay - dr1dot_dvy * Uyy ...
             + r2dot_yvz * az - dr1dot_dvz * Uzy;

% ∂(r1ddot)/∂z
dr1ddot_dz = r2dot_zx * vx + r2dot_zy * vy + r2dot_zz * vz ...
             + r2dot_zvx * ax - dr1dot_dvx * Uzx ...
             + r2dot_zvy * ay - dr1dot_dvy * Uzy ...
             + r2dot_zvz * az - dr1dot_dvz * Uzz;

% ∂(r1ddot)/∂vx
dr1ddot_dvx = r2dot_xvx * vx + dr1dot_dx ...
              + r2dot_yvx * vy + r2dot_zvx * vz ...
              - 2 * dr1dot_dvy;

% ∂(r1ddot)/∂vy
dr1ddot_dvy = r2dot_xvy * vx + dr1dot_dy ...
              + r2dot_yvy * vy + r2dot_zvy * vz ...
              + 2 * dr1dot_dvx;

% ∂(r1ddot)/∂vz
dr1ddot_dvz = r2dot_xvz * vx + r2dot_yvz * vy + r2dot_zvz * vz ...
              + dr1dot_dz;

end
