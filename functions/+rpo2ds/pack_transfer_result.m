function packed = pack_transfer_result(raw, target_name, target, settings)
%PACK_TRANSFER_RESULT Keep only the Run 02 data needed for figures and review.
%
% PROTOTYPE:
%   packed = rpo2ds.pack_transfer_result(raw, target_name, target, settings)
%
% DESCRIPTION:
%   Copy numeric trajectory, node, impulse, target, and solver summaries out
%   of the full transfer result. Missing optional payloads remain empty or
%   NaN; function handles and MATLAB objects are excluded from copied fields.
%
% INPUT:
%   raw          Full result from RPO2DS_transfer.
%   target_name  venus or mars.
%   target       Target struct with body_name and ECI/RTH v-infinity [km/s].
%   settings     Frozen validation settings to record with the compact result.
%
% OUTPUT:
%   packed       Struct with initial_guess, nodes, optimized trajectories,
%                optimized_impulses, metrics, convergence, and solver summaries.
%                Rotating coordinates/impulses use DU/VU; inertial positions
%                use km and DVvec_opt_in uses m/s. Metrics use m/s and days.
%                completed is based on a positive solver exitflag.

% Retain the target, patched seed, node layout, and optimized geometry.
packed = struct();
packed.target_name = string(target_name);
packed.target = struct( ...
    'body_name', target.body_name, ...
    'vinf_earth_eci_km_s', target.vinf_earth_eci_km_s(:), ...
    'vinf_earth_rth_km_s', target.vinf_earth_rth_km_s(:));
packed.initial_guess = local_initial_guess(raw);
packed.nodes = local_nodes(raw);
packed.optimized = local_trajectory(raw, 'trajectory');
packed.optimized_impulses = local_impulses(raw);
packed.metrics = struct( ...
    'total_delta_v_m_s', local_scalar(raw, 'dv_opt'), ...
    'time_of_flight_days', local_scalar(raw, 'tof_opt'), ...
    'earth_perigee_delta_v_m_s', local_earth_perigee_delta_v(raw), ...
    'rpo_departure_delta_v_m_s', local_rpo_departure_delta_v(raw), ...
    'terminal_vinf_eci_km_s', target.vinf_earth_eci_km_s(:), ...
    'terminal_vinf_rth_km_s', target.vinf_earth_rth_km_s(:));

% Keep solver diagnostics separately from physical result metrics.
packed.convergence = local_convergence(raw);
packed.solver_settings = settings;
packed.source_solver_fields = local_solver_summary(raw);
end

function initial_guess = local_initial_guess(raw)
%LOCAL_INITIAL_GUESS Copy patched arcs and an available preoptimization NLP trajectory.

initial_guess = struct();
if isfield(raw, 'initial_guess')
    source = raw.initial_guess;
    initial_guess = local_copy(source, {'xx', 'xxRot', 'xxFB', 'xxTot', 't', ...
        'tRot', 'tFB', 'tTot', 'xxF', 'xxB', 'xxBF'});
end
if isfield(raw, 'preopt') && isfield(raw.preopt, 'setup') && ~isempty(raw.preopt.setup)
    initial_guess.nlp = local_trajectory(raw, 'preopt');
end
end

function nodes = local_nodes(raw)
%LOCAL_NODES Copy node-layout metadata and initial shooting states when present.

nodes = struct();
if ~isfield(raw, 'nodes')
    return
end
source = raw.nodes;
nodes = local_copy(source, {'n_nodes', 'n_nodes_base', 'n_nodes_after_peri', ...
    'node_count_mode', 'nodes_per_day', 'transfer_time_days', ...
    'X_nodes_IG', 't_nodes_IG', 'xx_nodes_IG'});
end

function trajectory = local_trajectory(raw, field_name)
%LOCAL_TRAJECTORY Copy the retained preoptimization or optimized trajectory fields.

trajectory = struct();
if strcmp(field_name, 'preopt')
    source = raw.preopt.setup;
else
    if ~isfield(raw, field_name) || isempty(raw.(field_name))
        return
    end
    source = raw.(field_name);
end
trajectory = local_copy(source, {'t_nodes', 'xx_nodes', 'U_nodes', 'xx_traj_rot', ...
    't_traj', 'xx_traj_in', 'xx_nodes_in', 'DV_mag_opt', ...
    'DVvec_opt_rot', 'DVvec_opt_in', 'dv_opt', 'tof_opt', 'vInfMag'});
end

function value = local_scalar(source, field_name)
%LOCAL_SCALAR Return an available nonempty field or NaN for missing compact metrics.

if isfield(source, field_name) && ~isempty(source.(field_name))
    value = source.(field_name);
else
    value = NaN;
end
end

function value = local_earth_perigee_delta_v(raw)
%LOCAL_EARTH_PERIGEE_DELTA_V Read the final inertial impulse [m/s] or the patched perigee cost.

value = NaN;
if isfield(raw, 'dv') && isfield(raw.dv, 'DVvec_opt_in') && ...
        ~isempty(raw.dv.DVvec_opt_in)
    dv = vecnorm(raw.dv.DVvec_opt_in, 2, 2);
    value = dv(end);
    return
end
if isfield(raw, 'FullSolution') && istable(raw.FullSolution) && ...
        ismember('deltavpMag', raw.FullSolution.Properties.VariableNames)
    value = 1000 * raw.FullSolution.deltavpMag(1);
end
end

function value = local_rpo_departure_delta_v(raw)
%LOCAL_RPO_DEPARTURE_DELTA_V Read the initial inertial impulse [m/s] or scale the patched RPO burn.

value = NaN;
if isfield(raw, 'dv') && isfield(raw.dv, 'DVvec_opt_in') && ...
        ~isempty(raw.dv.DVvec_opt_in)
    dv = vecnorm(raw.dv.DVvec_opt_in, 2, 2);
    value = dv(1);
    return
end
if isfield(raw, 'FullSolution') && istable(raw.FullSolution) && ...
        ismember('deltaVRPO', raw.FullSolution.Properties.VariableNames)
    velocity_unit = NaN;
    if isfield(raw, 'constants') && isfield(raw.constants, 'VU')
        velocity_unit = raw.constants.VU;
    end
    if isfinite(velocity_unit)
        value = 1000 * velocity_unit * raw.FullSolution.deltaVRPO(1);
        return
    end
end
end

function impulses = local_impulses(raw)
%LOCAL_IMPULSES Retain normalized rotating impulses and inertial impulses [m/s].

impulses = struct();
if ~isfield(raw, 'dv') || isempty(raw.dv)
    return
end
if isfield(raw.dv, 'DV_mag_opt')
    impulses.magnitude_nd = raw.dv.DV_mag_opt(:);
end
if isfield(raw.dv, 'DVvec_opt_rot')
    impulses.rotating_nd = raw.dv.DVvec_opt_rot;
end
if isfield(raw.dv, 'DVvec_opt_in')
    impulses.inertial_m_s = raw.dv.DVvec_opt_in;
end
end

function convergence = local_convergence(raw)
%LOCAL_CONVERGENCE Collect solver diagnostics and mark positive exit status as completed.

convergence = struct('exitflag', NaN, 'objective', NaN, ...
    'iterations', NaN, 'constraint_residual_inf', NaN, 'completed', false);
if isfield(raw, 'exitflag') && ~isempty(raw.exitflag)
    convergence.exitflag = raw.exitflag;
end
if isfield(raw, 'J_opt') && ~isempty(raw.J_opt)
    convergence.objective = raw.J_opt;
end
if isfield(raw, 'output') && isstruct(raw.output)
    if isfield(raw.output, 'iterations')
        convergence.iterations = raw.output.iterations;
    end
    if isfield(raw.output, 'constrviolation')
        convergence.constraint_residual_inf = raw.output.constrviolation;
    end
end
convergence.completed = isfinite(convergence.exitflag) && convergence.exitflag > 0;
end

function summary = local_solver_summary(raw)
%LOCAL_SOLVER_SUMMARY Copy the active phase, geometry, timing, and target-frame labels.

summary = struct();
for field_name = {'operational_mode', 'target_frame', 'phase_mode', ...
        'time_mode', 'time_variable_mode', 'geometry_mode', ...
        'phase_seed', 'phase_final', 'seed_tof_exact_days', 'active_max_tof_days'}
    name = field_name{1};
    if isfield(raw, name)
        summary.(name) = raw.(name);
    end
end
end

function destination = local_copy(source, names)
%LOCAL_COPY Copy selected fields while excluding function handles and MATLAB objects.

destination = struct();
for index = 1:numel(names)
    name = names{index};
    if isfield(source, name)
        value = source.(name);
        if ~(isa(value, 'function_handle') || isobject(value))
            destination.(name) = value;
        end
    end
end
end
