function [J, dJdy] = objfun_CR3BP_regularized(y, problem)
%OBJFUN_CR3BP_REGULARIZED Minimize total normalized regularized impulse magnitude.
%
% PROTOTYPE:
%   [J, dJdy] = objfun_CR3BP_regularized(y, problem)
%
% DESCRIPTION:
%   For each active node, sum u^2 + w^2 + s^2. Under Oshima's quadratic impulse
%   mapping this equals norm(delta-v), not its square. The state/time/phase
%   entries do not contribute directly to the objective.
%
% INPUT:
%   y        Packed decision vector with first 9*N entries grouped as
%              [x, y, z, vx, vy, vz, u, w, s] per rotating-frame node.
%   problem  Struct with n and optional active_impulse_mask. Optional for the
%              legacy 10*N layout, where N is inferred and all burns are active.
%
% OUTPUT:
%   J        Scalar sum of normalized impulse magnitudes [VU].
%   dJdy     Optional [numel(y) x 1] analytic gradient; other entries are zero.

if nargin > 1 && isstruct(problem) && isfield(problem, 'n') && ~isempty(problem.n)
    N = problem.n;
else
    N = floor(length(y) / 10);
end
if nargin > 1
    active_impulse_mask = local_active_impulse_mask(problem, N);
else
    active_impulse_mask = true(N, 1);
end

% Initialize
J    = 0.0;
dJdy = [];

% Extract all u, w, s from y
% y layout:
%   node j -> y(9*(j-1)+1 : 9*(j-1)+9) = [x y z vx vy vz u w s]
%
u_all = zeros(N,1);
w_all = zeros(N,1);
s_all = zeros(N,1);

for j = 1:N
    base = 9*(j-1);
    u_all(j) = y(base + 7);
    w_all(j) = y(base + 8);
    s_all(j) = y(base + 9);
end

% Compute cost only at DV nodes
for j = 1:N
    if ~active_impulse_mask(j)
        continue
    end
    u = u_all(j);
    w = w_all(j);
    s = s_all(j);

    J = J + (u^2 + w^2 + s^2);
end

% If gradient is requested, compute analytic gradient
if nargout > 1
    dJdy = zeros(numel(y), 1);

    for j = 1:N
        if ~active_impulse_mask(j)
            continue
        end
        base = 9*(j-1);

        u = u_all(j);
        w = w_all(j);
        s = s_all(j);

        % d(u^2 + w^2 + s^2)/du = 2*u, etc.
        dJdy(base + 7) = dJdy(base + 7) + 2*u;
        dJdy(base + 8) = dJdy(base + 8) + 2*w;
        dJdy(base + 9) = dJdy(base + 9) + 2*s;
    end
end

end

function active_impulse_mask = local_active_impulse_mask(problem, n)
%LOCAL_ACTIVE_IMPULSE_MASK Read the active impulse mask and default to every node.

active_impulse_mask = true(n, 1);
if nargin > 0 && isstruct(problem) && isfield(problem, 'active_impulse_mask') && ~isempty(problem.active_impulse_mask)
    active_impulse_mask = logical(problem.active_impulse_mask(:));
end
end
