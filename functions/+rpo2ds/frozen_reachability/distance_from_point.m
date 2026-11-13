function [d_min, d_vec, is_below_min] = distance_from_point(states, point, d_min_req)
%DISTANCE_FROM_POINT  Distance of state(s) from a reference point.
%
%   [d_min, d_vec, is_below_min] = distance_from_point(states, point, d_min_req)
%
%   This function extracts positions from one or more state/position vectors
%   and computes the Euclidean distance to a given reference point.
%
%   Supported state formats (position-only or full state):
%     - 1x2, 1x3, 1x4, 1x5, 1x6   (single state/position, row)
%     - 2x1, 3x1, 4x1, 5x1, 6x1   (single state/position, column)
%     - Nx2, Nx3, Nx4–Nx6         (N states as rows)
%     - 2xN, 3xN, 4xN–6xN         (N states as columns)
%   Only the first 2 or 3 components are used as position:
%     - 2 components → [x, y, 0]
%     - ≥3 components → [x, y, z]
%
%   INPUT:
%       states    - State or position array in any of the formats above.
%
%       point     - Reference point. Accepted formats:
%                     * scalar: x_p       → [x_p, 0, 0]
%                     * 1x2 or 2x1: [x,y] → [x, y, 0]
%                     * 1x3 or 3x1: [x,y,z]
%
%       d_min_req - (optional) distance threshold. If provided, the third
%                   output is a boolean flag indicating whether
%                   d_min <= d_min_req. If not provided or empty, the
%                   third output is [].
%
%   OUTPUT:
%       d_min       - Minimum distance from 'point' over all states.
%
%       d_vec       - Column vector (N x 1) of distances, one per state.
%
%       is_below_min - Boolean flag:
%                        * true  if d_min <= d_min_req
%                        * false if d_min >  d_min_req
%                        * []    if d_min_req not provided/empty
%
    % 1) Normalize 'states' shape
    if isvector(states)
        % Single state: make it a 1xD row
        states = states(:).';  % force row
    end

    [nRows, nCols] = size(states);

    % Decide if states are stored by rows (N x D) or by columns (D x N)
    % We allow D in [2..6] on either dimension.
    if nCols >= 2 && nCols <= 6
        % States are N x D (rows)
        states_row = states;
    elseif nRows >= 2 && nRows <= 6
        % States are D x N (columns) -> transpose to N x D
        states_row = states.';  % N x D
        [nRows, nCols] = size(states_row);
    else
        error('distance_from_point:statesSize', ...
              'states must have 2–6 components in one dimension (rows or columns).');
    end

    % Use only first 2 or 3 components as position
    if nCols >= 3
        r = states_row(:, 1:3);             % [x y z]
    elseif nCols == 2
        r = [states_row, zeros(nRows, 1)];  % [x y 0]
    else
        error('distance_from_point:statesComponents', ...
              'states must have at least 2 components (x, y).');
    end

    N = size(r, 1);
    % 2) Normalize 'point'
    if isscalar(point)
        % Scalar: [x_p, 0, 0]
        r_point = [point, 0, 0];
    else
        p = point(:);  % column
        if numel(p) == 2
            r_point = [p(1), p(2), 0];      % [x, y, 0]
        elseif numel(p) == 3
            r_point = [p(1), p(2), p(3)];   % [x, y, z]
        else
            error('distance_from_point:pointSize', ...
                  'point must be scalar, 2-element, or 3-element vector.');
        end
    end
    % 3) Distance computation
    diff = r - repmat(r_point, N, 1);   % N x 3
    d_vec = sqrt(sum(diff.^2, 2));      % N x 1
    d_min = min(d_vec);
    % 4) Optional distance check
    if nargout > 2
        if nargin < 3 || isempty(d_min_req)
            is_below_min = [];
        else
            is_below_min = (d_min + 1e-5 <= d_min_req);
        end
    end

end

