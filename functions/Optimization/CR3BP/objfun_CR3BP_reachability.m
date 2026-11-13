function [J, dJdy] = objfun_CR3BP_reachability(y, problem)
%OBJFUN_CR3BP_REACHABILITY Maximize the terminal Earth-relative v-infinity magnitude.
%
% PROTOTYPE:
%   [J, dJdy] = objfun_CR3BP_reachability(y, problem)
%
% DESCRIPTION:
%   Minimize the negative final decision variable. The maneuver budget is
%   enforced by the nonlinear constraints, not by this objective.
%
% INPUT:
%   y        Packed reachability decision vector; y(end) is v-infinity [km/s].
%              The preceding state/time/phase entries follow the active solver layout.
%   problem  Optional node-count metadata; does not affect the objective.
%
% OUTPUT:
%   J        Negative reachable speed [km/s].
%   dJdy     Optional [numel(y) x 1] gradient, with -1 only at the last entry.

if nargin > 1 && isstruct(problem) && isfield(problem, 'n') && ~isempty(problem.n)
    N = problem.n;
else
    N = floor((length(y) - 1) / 10);
end

% Initialize
dJdy = [];

% Extract vInfMag from y
vInfMag = y(end);

J = -vInfMag;

% If gradient is requested, compute analytic gradient
if nargout > 1
    dJdy = zeros(numel(y), 1);

    dJdy(end) = -1;
end

end
