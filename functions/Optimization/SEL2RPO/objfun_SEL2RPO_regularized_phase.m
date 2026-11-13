function [J, dJdy] = objfun_SEL2RPO_regularized_phase(y, nlp)
%OBJFUN_SEL2RPO_REGULARIZED_PHASE Sum regularized impulse costs with their exact decision gradient.
%
% PROTOTYPE:
%   [J, dJdy] = objfun_SEL2RPO_regularized_phase(y, nlp)
%
% DESCRIPTION:
%   Compute J=sum_j(u_j^2+w_j^2+s_j^2). Under either supported
%   Oshima impulse orientation, this is the sum of impulse magnitudes
%   in normalized EM velocity units. State, duration and boundary-phase
%   variables do not contribute directly to the objective.
%
% INPUT:
%   y    [nvars x 1] Decision vector; frozen phase-boundary layout
%          has nvars=10*N-11: U1, interior (Xj,Uj), UN, segment
%          durations, SEL2 phase and RPO phase. Uj=[u,w,s] is mapped
%          through nlp.map.node(j).U; see the constraint function.
%   nlp  Struct with N (node count) and map.node(j).U (three indices).
%
% OUTPUT:
%   J     Scalar total delta-v [VU_EM]. Multiply by the problem
%           VU_km_s*1000 to express it in m/s.
%   dJdy  Optional [nvars x 1] gradient, zero outside impulse indices.

if nargin < 2 || ~isstruct(nlp)
    error('objfun_SEL2RPO_regularized_phase requires input struct nlp.');
end
if ~isfield(nlp, 'N') || ~isscalar(nlp.N)
    error('nlp.N is required and must be a scalar.');
end
if ~isfield(nlp, 'map') || ~isstruct(nlp.map) || ~isfield(nlp.map, 'node')
    error('nlp.map.node is required.');
end

% Read number of nodes and initialize outputs.
N = round(nlp.N);
J = 0.0;

% Preallocate gradient only when requested by caller.
if nargout > 1
    dJdy = zeros(numel(y), 1);
else
    dJdy = [];
end

% Accumulate cost and gradient node-by-node using mapped U indices.
for j = 1:N
    idxU = nlp.map.node(j).U;
    if isempty(idxU) || numel(idxU) ~= 3
        error('nlp.map.node(%d).U must contain 3 indices.', j);
    end

    u = y(idxU(1));
    w = y(idxU(2));
    s = y(idxU(3));

    J = J + (u^2 + w^2 + s^2);

    if nargout > 1
        dJdy(idxU(1)) = dJdy(idxU(1)) + 2*u;
        dJdy(idxU(2)) = dJdy(idxU(2)) + 2*w;
        dJdy(idxU(3)) = dJdy(idxU(3)) + 2*s;
    end
end
end
