function [Vrot, R] = rotation(V, axis, angle)
%ROTATION  Rotate 3D vector(s) about an axis by a given angle.
%
% USAGE
%   [Vrot, R] = rotation(V, axis, angle)
%
% INPUTS
%   V     : Vector or matrix of vectors to rotate.
%           Accepts: 1x3, 3x1, N×3, or 3×N. Each vector is [x y z].
%   axis  : Rotation axis vector (1x3 or 3x1). Will be normalized.
%   angle : Rotation angle in radians (scalar). Applied to all rows of V.
%
% OUTPUTS
%   Vrot  : Rotated vectors as an N×3 matrix (each row is one rotated vector).
%   R     : 3×3 rotation matrix used (Rodrigues' formula).
%
% NOTES
%   - If V is 3×N, it is interpreted as N column vectors and converted to N×3.
%   - If V is a single row/column vector, output is 3x1.
%   - If 'axis' is zero-norm and angle ≈ 0, no rotation is applied (R = I).
%     Otherwise a zero-norm axis is an error.
    axis = axis(:).';                 % make 1x3 row
    if numel(axis) ~= 3
        error('axis must have 3 elements.');
    end
    na = norm(axis);
    if na == 0
        if abs(angle) < eps
            R = eye(3);
        else
            error('Rotation axis must be non-zero unless angle is zero.');
        end
    else
        a = axis / na;                % unit axis
        % Skew-symmetric matrix of axis (Rodrigues' formula)
        K = [    0,   -a(3),  a(2);
               a(3),     0,  -a(1);
              -a(2),  a(1),     0 ];
        R = eye(3) + sin(angle)*K + (1 - cos(angle))*(K*K);
    end
    if isvector(V) && numel(V) == 3
        V = V(:).';                   % 1×3
    elseif size(V,2) == 3             % already N×3
        % keep as is
    elseif size(V,1) == 3             % 3×N -> transpose to N×3
        V = V.';
    else
        error('V must be 1x3, 3x1, N×3, or 3×N.');
    end
    % For row vectors, v_rot_row = v_row * R'
    Vrot = V * R.';
    if numel(Vrot) == 3
        Vrot = Vrot';
    end
end
