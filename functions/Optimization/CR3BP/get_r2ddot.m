function [r2, r2dot, r2ddot] = get_r2ddot(X, mu)
%GET_R2DDOT Evaluate distance and radial derivatives from the Moon.
%
% PROTOTYPE:
%   [r2, r2dot, r2ddot] = get_r2ddot(X, mu)
%
% DESCRIPTION:
%   Use the normalized 3D CR3BP rotating-frame acceleration to evaluate
%   r, dr/dt, and d^2r/dt^2 from the primary at x = 1 - mu. A positive
%   radial velocity means increasing distance from that primary.
%
% INPUT:
%   X   [6 elements] State [x, y, z, vx, vy, vz] in rotating DU/VU.
%   mu  Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   r2       Scalar primary distance [DU].
%   r2dot    Scalar radial velocity [DU/TU].
%   r2ddot   Scalar radial acceleration [DU/TU^2].

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

% Distance to the first primary (located at x = - mu)
dx1 = x - 1 + mu;

r2 = sqrt( dx1.^2 + y.^2 + z.^2 );

% First derivative: radial velocity
%
%   rdot = (r · v) / |r|
%
r2dot = (dx1 .* vx + y .* vy + z .* vz) ./ r2;

% Second derivative: radial acceleration
%
%   rddot = (v·v + r·a)/r  -  (r·v)*rdot / r²
%
%   Derived from:
%       r = |R|
%       rdot  = (R·V)/r
%       rddot = (V·V + R·A)/r  -  (R·V)*rdot / r²
%
term1 = (vx.^2 + vy.^2 + vz.^2 + dx1 .* ax + y .* ay + z .* az) ./ r2;
term2 = (dx1 .* vx + y .* vy + z .* vz) .* r2dot ./ (r2.^2);

r2ddot = term1 - term2;

end
