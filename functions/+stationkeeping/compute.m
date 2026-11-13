function result = compute(cfg,inputs)
%COMPUTE Run the two serial stationkeeping campaigns from frozen corrected states.
%
% PROTOTYPE:
%   result = stationkeeping.compute(cfg, inputs)
%
% DESCRIPTION:
%   Build both crossing maps before starting trials; no orbit solve
%   is performed. Reset the configured RNG at the start of each branch
%   and preserve the historical insertion, navigation and execution
%   error draw order. Keep failures and their diagnostics in the result;
%   compute cost statistics from finite completed trials only.
%
% INPUT:
%   cfg  Stationkeeping settings from the cfg.stationkeeping structure
%        returned by reproduction_config():
%        dynamics, units, reference, targeting, propagation, errors,
%        trial count, RNG settings, horizon and navigation delay.
%        See config/reproduction_config.m for the frozen field values.
%   inputs  Struct from data/inputs/run_05_reference_inputs.mat with
%             provenance and left/right branch structs. Each branch has
%             state_du [6 x 1] = [x;y;z;vx;vy;vz] in the rotating EM
%             barycentric frame (position [DU], velocity [VU]), and
%             sun_phase_rad [rad] at absolute time zero.
%
% OUTPUT:
%   result  Struct with settings, input_provenance, origin, created_at,
%             references.left/right (stationkeeping.reference) and
%             branches.left/right. Each campaign stores annual_cost_m_s,
%             completed and maneuver_count [K x 1], trials [K x 1]
%             (run_trial schema below), elapsed_s [s] and statistics
%             (stationkeeping.statistics). K is trials_per_branch.
%             Costs sum realized impulses before the configured horizon
%             [m/s]; for the frozen 365-day profile these are annual
%             costs, without rescaling for the final coast overshoot.

result = struct('settings',cfg,'input_provenance',inputs.provenance, ...
    'origin','Fresh public Run 05 computation','created_at',datestr(now,30));
branches = {'left','right'};
% Validate both sampled reference maps before starting either campaign.
for b = 1:numel(branches)
    name = branches{b};
    result.references.(name) = stationkeeping.reference(cfg,inputs.(name));
end
for b = 1:numel(branches)
    name = branches{b};
    reference = result.references.(name);

    % Restart the same random stream for each branch and allocate all trial records.
    rng(cfg.random_seed,cfg.random_generator);
    elapsed = tic;
    campaign = struct('annual_cost_m_s',nan(cfg.trials_per_branch,1), ...
        'completed',false(cfg.trials_per_branch,1), ...
        'maneuver_count',zeros(cfg.trials_per_branch,1));
    campaign.trials = repmat(empty_trial(),cfg.trials_per_branch,1);
    fprintf('\nStationkeeping %s: %d serial trials.\n',name,cfg.trials_per_branch);

    % Run serial trials without discarding incomplete cases or their diagnostics.
    for k = 1:cfg.trials_per_branch
        trial = run_trial(reference,cfg);
        campaign.trials(k) = trial;
        campaign.annual_cost_m_s(k) = trial.annual_cost_m_s;
        campaign.maneuver_count(k) = numel(trial.impulse_cost_m_s);
        campaign.completed(k) = strcmp(trial.stop_reason,'horizon_reached');
        if mod(k,25) == 0 || k == cfg.trials_per_branch || ~campaign.completed(k)
            fprintf('  %s %d/%d: %s; cost %.6f m/s; %d maneuvers; elapsed %.1f min\n', ...
                name,k,cfg.trials_per_branch,trial.stop_reason, ...
                campaign.annual_cost_m_s(k),campaign.maneuver_count(k),toc(elapsed)/60);
        end
    end

    % Summarize only finite costs from trials that reached the horizon.
    campaign.elapsed_s = toc(elapsed);
    campaign.statistics = stationkeeping.statistics(campaign.annual_cost_m_s,campaign.completed);
    result.branches.(name) = campaign;
    fprintf('%s finished: %d/%d complete; mean %.12f, P99 %.12f m/s; %.1f min.\n', ...
        name,sum(campaign.completed),cfg.trials_per_branch, ...
        campaign.statistics.mean_m_s,campaign.statistics.p99_m_s,campaign.elapsed_s/60);
end
end

function trial = run_trial(reference,cfg)
%RUN_TRIAL Simulate one noisy insertion and delayed-navigation correction sequence.
%
% PROTOTYPE:
%   trial = run_trial(reference, cfg)
%
% DESCRIPTION:
%   Coast the true state to a correction epoch. Reconstruct a
%   delayed observation, add navigation error, and propagate the estimate
%   forward before targeting the next reference crossing. Apply the
%   commanded impulse with independent componentwise execution error.
%   Carry coast-time corrections relative to the nominal reference
%   epoch. Stop on timing/targeting failure or after reaching the horizon;
%   a final coast can overshoot it, with no impulse applied afterward.
%
% INPUT:
%   reference  Crossing map from stationkeeping.reference.
%   cfg        Frozen stationkeeping configuration; see compute above.
%                All propagated states use [x;y;z;vx;vy;vz] in the EM
%                rotating barycentric frame (DU, VU); times use TU.
%
% OUTPUT:
%   trial  Struct containing annual_cost_m_s (realized cost sum [m/s]),
%           impulse_time_tu/impulse_cost_m_s [J x 1, TU and m/s],
%           target_time_tu [L x 1, correction epochs in TU],
%           requested_coast_tu/corrected_coast_tu [L x 1, TU],
%           target_success/target_iterations [L x 1],
%           target_stop_reason [L x 1 cell], and
%           target_residual_normalized [L x 4] in [x,y,z,vy] order.
%           Residuals describe the estimated commanded trajectory
%           before execution error, scaled by component tolerances.
%           stop_reason labels completion/failure; final_time_tu is the
%           last coast endpoint [TU]. J is applied impulses; L is
%           targeting attempts and can include one failed attempt.

trial = empty_trial();

% Draw the initial insertion error and start toward the first future crossing.
x_true = reference.initial_state_du+state_error(cfg);
t_now = 0;
[~,next_time] = point(reference,2);
coast_tu = next_time;
target_index = 1;
timing_offset_tu = 0;
year_tu = cfg.horizon_days/cfg.units.days_per_tu;
delay_tu = cfg.navigation_delay_days/cfg.units.days_per_tu;

% Advance the true state before deciding whether another correction is needed.
while true
    if coast_tu <= 0
        trial.stop_reason = 'non_monotonic_impulse_time';
        return;
    end
    x_end = coast(x_true,t_now,t_now+coast_tu,cfg,reference.sun_phase_rad);
    t_now = t_now+coast_tu;
    trial.final_time_tu = t_now;
    % The final coast can overshoot one year; never apply an impulse afterward.
    if t_now >= year_tu
        trial.stop_reason = 'horizon_reached';
        return;
    end

    % Model delayed navigation by back-propagating, perturbing, then repropagating.
    x_back = coast(x_end,t_now,t_now-delay_tu,cfg,reference.sun_phase_rad);
    x_estimated = x_back+state_error(cfg);
    estimate = coast(x_estimated,t_now-delay_tu,t_now,cfg,reference.sun_phase_rad);

    % Select the next nominal crossing and include the accumulated timing offset.
    target_index = target_index+1;
    [target_state,target_time] = point(reference,target_index+1);
    nominal_dt = target_time-t_now;
    requested_dt = nominal_dt+timing_offset_tu;
    if nominal_dt <= 0 || t_now+requested_dt <= t_now
        trial.stop_reason = 'non_monotonic_target_time';
        return;
    end
    % Preserve the historical subtraction used by the wrapper's dt bounds.
    target_end = t_now+requested_dt;
    solution = stationkeeping.target(estimate,t_now,target_end-t_now, ...
        target_state,cfg,reference.sun_phase_rad);

    % Retain every targeting attempt, including its last residual and stop reason.
    trial.target_time_tu(end+1,1) = t_now;
    trial.requested_coast_tu(end+1,1) = requested_dt;
    trial.corrected_coast_tu(end+1,1) = solution.coast_tu;
    trial.target_success(end+1,1) = solution.success;
    trial.target_stop_reason{end+1,1} = solution.stop_reason;
    trial.target_iterations(end+1,1) = solution.iterations;
    trial.target_residual_normalized(end+1,:) = solution.residual_normalized;
    if ~solution.success
        trial.stop_reason = ['targeting_' solution.stop_reason];
        return;
    end
    coast_tu = solution.coast_tu;
    if ~isfinite(coast_tu) || coast_tu <= 0
        trial.stop_reason = 'invalid_corrected_coast';
        return;
    end
    % Offset is relative to the nominal epoch, not to the previous request.
    timing_offset_tu = coast_tu-nominal_dt;

    % Apply execution error to the command and charge the realized impulse magnitude.
    execution_error = (cfg.errors.execution_sigma_fraction*ones(3,1)).*randn(3,1);
    dv_true = solution.dv_plan+solution.dv_plan.*execution_error;
    dv_m_s = norm(dv_true,2)*cfg.units.VU_km_s*1000;
    x_end(4:6) = x_end(4:6)+dv_true;
    trial.impulse_time_tu(end+1,1) = t_now;
    trial.impulse_cost_m_s(end+1,1) = dv_m_s;
    trial.annual_cost_m_s = trial.annual_cost_m_s+dv_m_s;
    x_true = x_end;
end
end

function trial = empty_trial()
%EMPTY_TRIAL Allocate an empty trial record; field meanings and units are in run_trial.

trial = struct('annual_cost_m_s',0,'impulse_time_tu',zeros(0,1),'impulse_cost_m_s',zeros(0,1), ...
    'target_time_tu',zeros(0,1),'requested_coast_tu',zeros(0,1), ...
    'corrected_coast_tu',zeros(0,1),'target_success',false(0,1), ...
    'target_stop_reason',{cell(0,1)},'target_iterations',zeros(0,1), ...
    'target_residual_normalized',zeros(0,4),'stop_reason','', 'final_time_tu',0);
end

function [state,time] = point(reference,index)
%POINT Return the initial anchor once, then repeat two future crossings per cycle.
% index is a positive sequence index. state is [6 x 1] in EM DU/VU;
% time is absolute TU. reference is the stationkeeping.reference map.

if index == 1
    state = reference.initial_state_du;
    time = 0;
else
    periodic_index = index-2;
    slot = mod(periodic_index,2)+1;
    time = floor(periodic_index/2)*reference.cycle_tu+reference.crossing_times_tu(slot);
    state = reference.crossing_states_du(slot,:).';
end
end

function error_state = state_error(cfg)
%STATE_ERROR Draw independent Gaussian position/velocity errors in preserved RNG order.
% Return [6 x 1] normalized EM errors: position sigma [km]/DU_km and
% velocity sigma [cm/s]*1e-5/VU_km_s, three components of each.

position = ((cfg.errors.position_sigma_km/cfg.units.DU_km)*ones(3,1)).*randn(3,1);
velocity = (((cfg.errors.velocity_sigma_cm_s*1e-5)/cfg.units.VU_km_s)*ones(3,1)).*randn(3,1);
error_state = [position;velocity];
end

function state = coast(x0,t0,tf,cfg,sun_phase)
%COAST Propagate a [6 x 1] EM rotating state from t0 to tf [TU], forward or backward.
% Use cfg dynamics.mu and propagation tolerances; sun_phase [rad] is
% defined at absolute time zero. Return state in [x;y;z;vx;vy;vz], DU/VU.

out = propagate(struct('x0',x0,'t0',t0,'tf',tf,'mu',cfg.dynamics.mu, ...
    'dimensions','3D','framework','BCR4BP','stm_computation','false', ...
    'thetaSun0',sun_phase,'absTol',cfg.propagation.absolute_tolerance, ...
    'relTol',cfg.propagation.relative_tolerance));
state = out.xf;
end
