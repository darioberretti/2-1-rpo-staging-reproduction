function [J, dJdy] = objfun_CR3BP_regularized_phase(y, problem)
%OBJFUN_CR3BP_REGULARIZED_PHASE Minimize total normalized regularized impulse magnitude.
%
% PROTOTYPE:
%   [J, dJdy] = objfun_CR3BP_regularized_phase(y, problem)
%
% DESCRIPTION:
%   For each active node, sum u^2 + w^2 + s^2. Under Oshima's quadratic impulse
%   mapping this equals norm(delta-v), not its square. The state/time/phase
%   entries do not contribute directly to the objective.
%
% INPUT:
%   y        Packed decision vector with first 9*N entries grouped as
%              [x, y, z, vx, vy, vz, u, w, s] per rotating-frame node.
%   problem  Struct with n and optional active_impulse_mask. Required for
%              the phase-aware layout.
%
% OUTPUT:
%   J        Scalar sum of normalized impulse magnitudes [VU].
%   dJdy     Optional [numel(y) x 1] analytic gradient; other entries are zero.

n = problem.n;
active_impulse_mask = local_active_impulse_mask(problem, n);
J = 0.0;
dJdy = [];

% Accumulate the regularized delta-v cost on every node.
for j = 1:n
    if ~active_impulse_mask(j)
        continue
    end
    base = 9 * (j - 1);
    u = y(base + 7);
    w = y(base + 8);
    s = y(base + 9);
    J = J + (u^2 + w^2 + s^2);
end

% Build the analytic gradient when the caller requests it.
if nargout > 1
    dJdy = zeros(numel(y), 1);
    for j = 1:n
        if ~active_impulse_mask(j)
            continue
        end
        base = 9 * (j - 1);
        dJdy(base + 7) = 2 * y(base + 7);
        dJdy(base + 8) = 2 * y(base + 8);
        dJdy(base + 9) = 2 * y(base + 9);
    end
end
end

function active_impulse_mask = local_active_impulse_mask(problem, n)
%LOCAL_ACTIVE_IMPULSE_MASK Read the active impulse mask and default to every node.

active_impulse_mask = true(n, 1);
if isfield(problem, 'active_impulse_mask') && ~isempty(problem.active_impulse_mask)
    active_impulse_mask = logical(problem.active_impulse_mask(:));
end
end
