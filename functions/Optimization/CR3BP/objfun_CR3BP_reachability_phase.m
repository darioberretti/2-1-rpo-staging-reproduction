function [J, dJdy] = objfun_CR3BP_reachability_phase(y, ~)
%OBJFUN_CR3BP_REACHABILITY_PHASE Maximize the terminal Earth-relative v-infinity magnitude.
%
% PROTOTYPE:
%   [J, dJdy] = objfun_CR3BP_reachability_phase(y, problem)
%
% DESCRIPTION:
%   Minimize the negative final decision variable. The maneuver budget is
%   enforced by the nonlinear constraints, not by this objective.
%
% INPUT:
%   y        Packed reachability decision vector; y(end) is v-infinity [km/s].
%              The preceding state/time/phase entries follow the active solver layout.
%   problem  Unused positional argument retained by the shared interface.
%
% OUTPUT:
%   J        Negative reachable speed [km/s].
%   dJdy     Optional [numel(y) x 1] gradient, with -1 only at the last entry.

J = -y(end);
dJdy = [];

% Differentiate only with respect to the reachable vInf magnitude.
if nargout > 1
    dJdy = zeros(numel(y), 1);
    dJdy(end) = -1;
end
end
