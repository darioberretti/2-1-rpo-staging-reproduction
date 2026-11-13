function [J, dJdy] = objfun_vleo_rpo_regularized_2d_fast(y, nlp)
%OBJFUN_VLEO_RPO_REGULARIZED_2D_FAST Sum active planar regularized impulse magnitudes.
%
% PROTOTYPE:
%   [J, dJdy] = objfun_vleo_rpo_regularized_2d_fast(y, nlp)
%
% DESCRIPTION:
%   Accumulate u^2+w^2 at the active nodes, which equals normalized impulse
%   magnitude under the planar Oshima map. This uses the same cost as the
%   legacy backend; acceleration is in the corresponding constraint pipeline.
%
% INPUT:
%   y    Mapped decision-vector column from the multiple-shooting setup.
%   nlp  setup.nlp; N, map.node(j).U, and active_impulse_mask select terms.
%
% OUTPUT:
%   J     Total normalized impulse magnitude [VU].
%   dJdy  Optional [numel(y) x 1] analytic gradient. State/time/phase and
%           inactive impulse entries have zero direct cost derivative.

N = nlp.N;
J = 0;
dJdy = [];

% Preallocate the gradient only when the caller requests it.
if nargout > 1
    dJdy = zeros(numel(y), 1);
end

% Accumulate the active node-wise impulse magnitudes.
for j = 1:N
    % Skip inactive impulse nodes in endpoint-only mode.
    if ~nlp.active_impulse_mask(j)
        continue
    end

    % Read the regularized variables at node j.
    idx_u = nlp.map.node(j).U;
    u = y(idx_u(1));
    w = y(idx_u(2));

    % Add the planar Oshima impulse magnitude.
    J = J + (u^2 + w^2);

    % Differentiate the node contribution analytically.
    if nargout > 1
        dJdy(idx_u(1)) = dJdy(idx_u(1)) + 2 * u;
        dJdy(idx_u(2)) = dJdy(idx_u(2)) + 2 * w;
    end
end

end
