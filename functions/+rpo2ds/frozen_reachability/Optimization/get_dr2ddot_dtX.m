function [dr2dot_dt, dr2dot_dx, dr2dot_dy, dr2dot_dz, ...
          dr2dot_dvx, dr2dot_dvy, dr2dot_dvz, ...
          dr2ddot_dt, dr2ddot_dx, dr2ddot_dy, dr2ddot_dz, ...
          dr2ddot_dvx, dr2ddot_dvy, dr2ddot_dvz] = get_dr2ddot_dtX(X, constants)
%GET_DR2DDOT_DTX Differentiate primary-2 radial velocity and acceleration.
%
% PROTOTYPE:
%   [dr2dot_dt, dr2dot_dx, dr2dot_dy, dr2dot_dz, ...
%       dr2dot_dvx, dr2dot_dvy, dr2dot_dvz, ...
%       dr2ddot_dt, dr2ddot_dx, dr2ddot_dy, dr2ddot_dz, ...
%       dr2ddot_dvx, dr2ddot_dvy, dr2ddot_dvz] = get_dr2ddot_dtX(X, constants)
%
% DESCRIPTION:
%   Evaluate analytic partials of dr2/dt and d^2r2/dt^2 with respect
%   to explicit time and the six normalized rotating-frame state components.
%   Explicit time derivatives are zero for the autonomous CR3BP.
%
% INPUT:
%   X          [6 elements] State [x, y, z, vx, vy, vz] in DU/VU.
%   constants  Legacy vector whose first entry is the mass ratio mu [-].
%
% OUTPUT:
%   First 7 scalars   Partials of r2dot with respect to
%                       [t, x, y, z, vx, vy, vz], in that order.
%   Last 7 scalars    Corresponding partials of r2ddot. Each derivative
%                       uses the output/input units from get_r2ddot.

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

% Distance factors for primary 2 (the one we care about for r2, r2dot, r2ddot)
r2_05 = (dx2^2 + y^2 + z^2)^(-0.5);   % 1 / |r2|
r2_15 = (dx2^2 + y^2 + z^2)^(-1.5);   % 1 / |r2|^3
r2_25 = (dx2^2 + y^2 + z^2)^(-2.5);   % 1 / |r2|^5

% Distance factors for primary 1 (still needed for second derivatives of U)
r1_15 = (dx1^2 + y^2 + z^2)^(-1.5);   % 1 / |r1|^3
r1_25 = (dx1^2 + y^2 + z^2)^(-2.5);   % 1 / |r1|^5

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

% First derivatives of r2dot

% r2dot does not explicitly depend on time
dr2dot_dt  = 0;

% ∂(r2dot)/∂x, ∂(r2dot)/∂y, ∂(r2dot)/∂z
dr2dot_dx = r2_15 * ( vx*(y^2 + z^2)   - dx2*(y*vy + z*vz) );
dr2dot_dy = r2_15 * ( vy*(dx2^2 + z^2) - y*(dx2*vx + z*vz) );
dr2dot_dz = r2_15 * ( vz*(dx2^2 + y^2) - z*(dx2*vx + y*vy) );

% ∂(r2dot)/∂vx, ∂(r2dot)/∂vy, ∂(r2dot)/∂vz
dr2dot_dvx = dx2 * r2_05;
dr2dot_dvy = y   * r2_05;
dr2dot_dvz = z   * r2_05;

% Second derivatives needed for r2ddot

% r2dot_xx, r2dot_xy, r2dot_xz, ...
r2dot_xx = -3 * r2_25 * dx2 * ( vx*(y^2 + z^2) - dx2*(y*vy + z*vz) ) ...
           - r2_15 * (y*vy + z*vz);

r2dot_xy = -3 * r2_25 * y * ( vx*(y^2 + z^2) - dx2*(y*vy + z*vz) ) ...
           + r2_15 * (2*y*vx - dx2*vy);

r2dot_xz = -3 * r2_25 * z * ( vx*(y^2 + z^2) - dx2*(y*vy + z*vz) ) ...
           + r2_15 * (2*z*vx - dx2*vz);

r2dot_yx = r2dot_xy;

r2dot_yy = -3 * r2_25 * y * ( vy*(dx2^2 + z^2) - y*(dx2*vx + z*vz) ) ...
           - r2_15 * (dx2*vx + z*vz);

r2dot_yz = -3 * r2_25 * z * ( vy*(dx2^2 + z^2) - y*(dx2*vx + z*vz) ) ...
           + r2_15 * (2*z*vy - y*vz);

r2dot_zx = r2dot_xz;
r2dot_zy = r2dot_yz;

r2dot_zz = -3 * r2_25 * z * ( vz*(dx2^2 + y^2) - z*(dx2*vx + y*vy) ) ...
           - r2_15 * (dx2*vx + y*vy);

% Mixed derivatives with velocities
r2dot_xvx =  r2_15 * (y^2 + z^2);
r2dot_xvy = -r2_15 * dx2 * y;
r2dot_xvz = -r2_15 * dx2 * z;

r2dot_yvx = -r2_15 * dx2 * y;
r2dot_yvy =  r2_15 * (dx2^2 + z^2);
r2dot_yvz = -r2_15 * y * z;

r2dot_zvx = -r2_15 * dx2 * z;
r2dot_zvy = -r2_15 * y * z;
r2dot_zvz =  r2_15 * (dx2^2 + y^2);

% Derivatives of r2ddot

% r2ddot does not explicitly depend on time in this formulation
dr2ddot_dt = 0;

% ∂(r2ddot)/∂x
dr2ddot_dx = r2dot_xx * vx + r2dot_xy * vy + r2dot_xz * vz ...
             + r2dot_xvx * ax - dr2dot_dvx * Uxx ...
             + r2dot_xvy * ay - dr2dot_dvy * Uyx ...
             + r2dot_xvz * az - dr2dot_dvz * Uzx;

% ∂(r2ddot)/∂y
dr2ddot_dy = r2dot_yx * vx + r2dot_yy * vy + r2dot_yz * vz ...
             + r2dot_yvx * ax - dr2dot_dvx * Uyx ...
             + r2dot_yvy * ay - dr2dot_dvy * Uyy ...
             + r2dot_yvz * az - dr2dot_dvz * Uzy;

% ∂(r2ddot)/∂z
dr2ddot_dz = r2dot_zx * vx + r2dot_zy * vy + r2dot_zz * vz ...
             + r2dot_zvx * ax - dr2dot_dvx * Uzx ...
             + r2dot_zvy * ay - dr2dot_dvy * Uzy ...
             + r2dot_zvz * az - dr2dot_dvz * Uzz;

% ∂(r2ddot)/∂vx
dr2ddot_dvx = r2dot_xvx * vx + dr2dot_dx ...
              + r2dot_yvx * vy + r2dot_zvx * vz ...
              - 2 * dr2dot_dvy;

% ∂(r2ddot)/∂vy
dr2ddot_dvy = r2dot_xvy * vx + dr2dot_dy ...
              + r2dot_yvy * vy + r2dot_zvy * vz ...
              + 2 * dr2dot_dvx;

% ∂(r2ddot)/∂vz
dr2ddot_dvz = r2dot_xvz * vx + r2dot_yvz * vy + r2dot_zvz * vz ...
              + dr2dot_dz;

end
