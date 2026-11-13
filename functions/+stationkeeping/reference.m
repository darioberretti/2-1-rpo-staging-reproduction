function reference = reference(cfg,input)
%REFERENCE Build a two-crossing reference map from the frozen branch state.
%
% PROTOTYPE:
%   reference = stationkeeping.reference(cfg, input)
%
% DESCRIPTION:
%   Propagate one configured two-period cycle on the historical
%   sample grid. Combine near-axis samples with linearly interpolated
%   sign changes, merge neighboring epochs, and retain the first two
%   future y=0, vy>0 crossings, including the cycle endpoint. This
%   sampled map is not an event-refined or newly corrected orbit.
%
% INPUT:
%   cfg  Stationkeeping settings from the cfg.stationkeeping structure
%        returned by reproduction_config():
%        dynamics, units, reference, targeting, propagation, errors,
%        trial count, RNG settings, horizon and navigation delay.
%        See config/reproduction_config.m for the frozen field values.
%   input  Frozen branch struct with state_du [6 x 1] and sun_phase_rad
%            [rad] at absolute time zero. State order is [x;y;z;vx;vy;vz]
%            in the EM rotating barycentric frame (position DU, velocity VU).
%
% OUTPUT:
%   reference  Struct with initial_state_du [6 x 1], sun_phase_rad,
%                cycle_tu [TU], crossing_times_tu [2 x 1, TU] and
%                crossing_states_du [2 x 6, EM DU/VU]. Crossing times
%                are relative to the cycle start. Fewer than two
%                admissible crossings raises an error before trials run.

% Sample the frozen BCR4BP orbit without terminal events or orbit correction.
out = propagate(struct('x0',input.state_du,'t0',0,'tf',cfg.reference.cycle_tu, ...
    'mu',cfg.dynamics.mu,'dimensions','3D','framework','BCR4BP', ...
    'stm_computation','false','thetaSun0',input.sun_phase_rad, ...
    'absTol',cfg.propagation.absolute_tolerance, ...
    'relTol',cfg.propagation.relative_tolerance, ...
    'tspan',linspace(0,cfg.reference.cycle_tu,cfg.reference.sample_count).'));

% Collect near-axis samples, then interpolate every sampled sign change.
tt = out.tt(:);
xx = out.xx(:,1:6);
on_axis = abs(xx(:,2)) < cfg.reference.y_crossing_tolerance;
crossing_times = tt(on_axis);
crossing_states = xx(on_axis,:);
for k = 1:numel(tt)-1
    y1 = xx(k,2); y2 = xx(k+1,2);
    if y1*y2 < 0
        fraction = y1/(y1-y2);
        crossing_times(end+1,1) = tt(k)+fraction*(tt(k+1)-tt(k)); %#ok<AGROW>
        crossing_states(end+1,:) = xx(k,:)+fraction*(xx(k+1,:)-xx(k,:)); %#ok<AGROW>
    end
end

% Order and merge candidate crossings before selecting future positive-vy ones.
[crossing_times,order] = sort(crossing_times);
crossing_states = crossing_states(order,:);
merge_tolerance = cfg.reference.time_merge_tolerance_tu;
keep = [true; abs(diff(crossing_times)) >= merge_tolerance];
crossing_times = crossing_times(keep);
crossing_states = crossing_states(keep,:);
keep = crossing_times > merge_tolerance & ...
    crossing_times <= cfg.reference.cycle_tu+merge_tolerance & crossing_states(:,5) > 0;
crossing_times = crossing_times(keep);
crossing_states = crossing_states(keep,:);
crossing_times(crossing_times > cfg.reference.cycle_tu) = cfg.reference.cycle_tu;
% The map must be usable before beginning either expensive trial campaign.
if numel(crossing_times) < 2
    error('stationkeeping:Crossings','The frozen reference needs two future positive-vy crossings.');
end

% Keep the initial anchor separate from the two periodically repeated crossings.
reference = struct('initial_state_du',xx(1,:).','sun_phase_rad',input.sun_phase_rad, ...
    'cycle_tu',cfg.reference.cycle_tu,'crossing_times_tu',crossing_times(1:2), ...
    'crossing_states_du',crossing_states(1:2,:));
end
