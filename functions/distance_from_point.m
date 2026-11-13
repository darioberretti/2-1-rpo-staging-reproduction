function [minimum_distance, distances, is_below_minimum] = ...
        distance_from_point(states, point, requested_minimum)
%DISTANCE_FROM_POINT Measure Cartesian state-row distances from a fixed point.
%
% PROTOTYPE:
%   [minimum_distance, distances, is_below_minimum] = distance_from_point( ...
%       states, point, requested_minimum)
%
% DESCRIPTION:
%   Use the first three state columns as positions. Accept a scalar point
%   as an x-axis location. The optional threshold flag uses the implemented
%   absolute 1e-5 distance tolerance; no unit conversion is applied.
%
% INPUT:
%   states             State vector or [N x M] numeric array with M >= 3.
%   point              [3 elements] Position, or scalar x coordinate.
%   requested_minimum  Optional distance threshold in the same position units.
%
% OUTPUT:
%   minimum_distance   Minimum sampled distance.
%   distances          [N x 1] Distance per row.
%   is_below_minimum   minimum_distance + 1e-5 <= requested_minimum,
%                        or [] when the threshold is absent/empty.

if isvector(states)
    states = states(:).';
end
if size(states, 2) < 3
    error('distance_from_point:StatesSize', ...
        'states must contain at least three position components per row.');
end
if isscalar(point)
    point = [point, 0, 0];
else
    point = point(:).';
end
if numel(point) ~= 3
    error('distance_from_point:PointSize', ...
        'point must be scalar or a three-component vector.');
end
distances = vecnorm(states(:, 1:3) - point, 2, 2);
minimum_distance = min(distances);
if nargout > 2
    if nargin < 3 || isempty(requested_minimum)
        is_below_minimum = [];
    else
        is_below_minimum = minimum_distance + 1e-5 <= requested_minimum;
    end
end
end
