function [f, gradf] = obj_zero(x)
%OBJ_ZERO Constant objective function for fmincon.
%
% DESCRIPTION
%   This function can be used as an objective function in fmincon when
%   you only care about satisfying constraints. It always returns f = 0.
%
% INPUT
%   x      - Decision variable vector (any size).
%
% OUTPUT
%   f      - Objective function value (always 0).
%   gradf  - Gradient of the objective function with respect to x.
%            This is a zero vector with the same size as x.
%
% NOTE
%   To use the gradient in fmincon, set:
%   options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true);

    % Constant objective value
    f = 0;

    % Return gradient if requested
    if nargout > 1
        gradf = zeros(size(x));
    end
end
