function [J, dJdy] = objfun_CR3BP_reachability(y)
%OBJFUN_CR3BP_REACHABILITY Maximize the terminal Earth-relative v-infinity magnitude.
%
% PROTOTYPE:
%   [J, dJdy] = objfun_CR3BP_reachability(y)
%
% DESCRIPTION:
%   Minimize the negative final decision variable. The maneuver budget is
%   enforced by the nonlinear constraints, not by this objective.
%
% INPUT:
%   y        Packed reachability decision vector; y(end) is v-infinity [km/s].
%              The historical layout is [9*N node entries; N times; vInfMag].
%
% OUTPUT:
%   J        Negative reachable speed [km/s].
%   dJdy     Optional [numel(y) x 1] gradient, with -1 only at the last entry.

nn = 10; % Number of variables per node
N  = (length(y) - 1)/nn; % Number of nodes

% Initialize
dJdy = [];

% Extract vInfMag from y
vInfMag = y(end);

J = -vInfMag;

% If gradient is requested, compute analytic gradient
if nargout > 1
    dJdy = zeros(10*N + 1, 1);

    dJdy(end) = -1;
end

end
