function [states, times, valid] = propagate_crossings(initial_state, count, time_limit, mu, compute_stm, event_tolerance, event_direction, reset_time)
%PROPAGATE_CROSSINGS Propagate through successive y = 0 crossings.
%
% PROTOTYPE:
%   [states, times, valid] = rpo.propagate_crossings(initial_state, count, ...
%       time_limit, mu, compute_stm, event_tolerance, event_direction, reset_time)
%
% DESCRIPTION:
%   Propagate planar CR3BP segments in the normalized Earth-Moon rotating
%   frame. Suppress IntegrationTolNotMet warnings only during this call.
%   Each requested segment must terminate at an accepted crossing.
%
% INPUT:
%   initial_state    [4 elements] State [x, y, vx, vy] in DU and VU.
%   count            Number of successive crossings to retain.
%   time_limit       Integration end time [TU] for each propagate call.
%   mu               Secondary-to-total mass ratio [-].
%   compute_stm      Logical; integrate a segment STM when true.
%   event_tolerance  Maximum accepted absolute y at a crossing [DU].
%   event_direction  Optional crossing direction: -1, 0, or 1; default 0.
%   reset_time       Optional logical, default false. If true, restart each
%                      segment at t = 0; otherwise use the previous end time.
%
% OUTPUT:
%   states           [count x 1] Cell array of state / stm structs. Each STM
%                      maps perturbations from its segment start to its end.
%   times            [count x 1] Crossing times [TU], accumulated when
%                      reset_time is true. Unfilled entries remain NaN.
%   valid            True only when every requested crossing is accepted.

% Restore the caller's warning state on every return path.
warning_state = warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
restore_warning = onCleanup(@() warning(warning_state));
if nargin < 7
    event_direction = 0;
end
if nargin < 8
    reset_time = false;
end

% Initialize the crossing chain and the requested STM mode.
states = cell(count, 1);
times = NaN(count, 1);
valid = false;
state = initial_state(:);
start_time = 0;
elapsed_time = 0;
if compute_stm
    stm_flag = 'true';
else
    stm_flag = 'false';
end

% Integrate successive segments and reject incomplete or off-section endpoints.
for crossing = 1:count
    [final_state, stm, ~, ~, ~, final_time] = propagate( ...
        state, start_time, time_limit, mu, '2D', 'CR3BP', stm_flag, ...
        true, event_direction);
    if final_time >= time_limit || abs(final_state(2)) > event_tolerance
        return
    end
    states{crossing} = struct('state', final_state, 'stm', stm);
    state = final_state;

    % Keep cumulative crossing times when restarting the integration clock.
    if reset_time
        elapsed_time = elapsed_time + final_time;
        times(crossing) = elapsed_time;
        start_time = 0;
    else
        times(crossing) = final_time;
        start_time = final_time;
    end
end
valid = true;
end
