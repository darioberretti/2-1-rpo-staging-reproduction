function dec = decode_vleo_rpo_decision_vector_2d(y, nlp)
%DECODE_VLEO_RPO_DECISION_VECTOR_2D Decode the mapped planar transfer decision vector.
%
% PROTOTYPE:
%   dec = decode_vleo_rpo_decision_vector_2d(y, nlp)
%
% DESCRIPTION:
%   Use nlp.map rather than assuming uniform state blocks: every node owns
%   two regularized impulse variables, while only interior nodes own four
%   state variables. Reconstruct absolute node epochs from segment durations.
%
% INPUT:
%   y    [nvars x 1] Decision vector from the multiple-shooting setup.
%   nlp  setup.nlp with N, map, and the fixed departure-altitude setting.
%
% OUTPUT:
%   dec  Struct with U_nodes [N x 2], X_nodes_interior [N x 4] (endpoint rows
%          remain NaN), dt_nodes [N-1 x 1], and t_nodes [N x 1] starting at 0.
%          States use rotating [x, y, vx, vy] in DU/VU; times and phase_rpo
%          use TU, phase_departure/theta_sun0 rad, departure_altitude_km km.

N = nlp.N;

U_nodes = zeros(N, 2);

X_nodes_interior = NaN(N, 4);

% Decode all node-wise regularized variables.
for j = 1:N
    idx_u = nlp.map.node(j).U;
    U_nodes(j, :) = y(idx_u(:)).';

    % Decode the interior state only when that node owns one.
    idx_x = nlp.map.node(j).X;
    if ~isempty(idx_x)
        X_nodes_interior(j, :) = y(idx_x(:)).';
    end
end

% Decode the segment durations.
dt_nodes = y(nlp.map.dt_nodes(:));

% Decode the scalar phase variables.
phase_departure = y(nlp.map.phase_departure);
phase_rpo = y(nlp.map.phase_rpo);
theta_sun0 = y(nlp.map.theta_sun0);

% Decode the departure altitude according to the selected mode.
if isempty(nlp.map.departure_altitude_km)
    % Use the fixed departure altitude when it is not a decision variable.
    departure_altitude_km = nlp.departure_altitude_km_fixed;
else
    % Read the departure altitude from the decision vector.
    departure_altitude_km = y(nlp.map.departure_altitude_km);
end

% Package the decoded decision variables.
dec = struct;
dec.U_nodes = U_nodes;
dec.X_nodes_interior = X_nodes_interior;
dec.dt_nodes = dt_nodes;
dec.phase_departure = phase_departure;
dec.phase_rpo = phase_rpo;
dec.theta_sun0 = theta_sun0;
dec.departure_altitude_km = departure_altitude_km;
dec.t_nodes = [0; cumsum(dt_nodes(:))];

end
