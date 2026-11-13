function C = jacobi_constant(state, mu)
%JACOBI_CONSTANT CR3BP Jacobi constant for planar or spatial states.
%
% PROTOTYPE:
%   C = rpo.jacobi_constant(state, mu)
%
% DESCRIPTION:
%   Evaluate C = 2U - v^2 in the normalized rotating barycentric frame,
%   with primaries at x = -mu and x = 1 - mu.
%
% INPUT:
%   state   [N x 4] Rows [x, y, vx, vy], or [N x 6] rows
%             [x, y, z, vx, vy, vz]. Positions [DU], velocities [VU].
%   mu      Secondary-to-total mass ratio [-].
%
% OUTPUT:
%   C       [N x 1] Jacobi constants [-].

% Interpret planar states as spatial states with zero out-of-plane components.
x = state(:, 1);
y = state(:, 2);
if size(state, 2) == 4
    z = zeros(size(x));
    velocity = state(:, 3:4);
else
    z = state(:, 3);
    velocity = state(:, 4:6);
end

% Evaluate the effective potential and subtract the squared rotating-frame speed.
r1 = sqrt((x + mu).^2 + y.^2 + z.^2);
r2 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);
C = x.^2 + y.^2 + 2 * (1 - mu) ./ r1 + 2 * mu ./ r2 ...
    - sum(velocity.^2, 2);
end
