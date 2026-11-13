function [f, gradf] = obj_zero(x)
%OBJ_ZERO Return a constant zero objective and optional zero gradient.
%
% PROTOTYPE:
%   [f, gradf] = obj_zero(x)
%
% DESCRIPTION:
%   Use for fmincon feasibility stages with constraints supplied
%   separately. Enable SpecifyObjectiveGradient to request gradf.
%
% INPUT:
%   x      Decision array; normally the fmincon decision column vector.
%
% OUTPUT:
%   f      Scalar zero.
%   gradf  Optional zero array with exactly the same size as x.

    f = 0;

    % Return gradient if requested
    if nargout > 1
        gradf = zeros(size(x));
    end
end
