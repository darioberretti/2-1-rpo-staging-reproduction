function solve_out = solve_casadi_ipopt_two_stage(problem_definition, y0, opts)
%SOLVE_CASADI_IPOPT_TWO_STAGE Run the IPOPT feasibility and objective stages.
%
% PROTOTYPE:
%   solve_out = solve_casadi_ipopt_two_stage(problem_definition, y0, opts)
%
% DESCRIPTION:
%   Optionally solve the zero-objective feasibility problem, then warm-start
%   the configured objective solve from that result. Map IPOPT return status
%   into exit flags while retaining raw statistics for both stages.
%
% INPUT:
%   problem_definition  CasADi NLP, variable/constraint bounds from
%                         build_casadi_rpo2ds_initial_state_problem.
%   y0                  Initial decision-vector column in the same layout.
%   opts                Iteration/tolerance/display settings, optional
%                         casadi_path, run_stage1 (default true), progress_logs.
%
% OUTPUT:
%   solve_out           Struct with y_stage1/y_opt, J_opt, final exitflag/output,
%                         stage-specific flags/statistics, and solver options.
%                         Physical decoding and independent residual checks
%                         belong to the caller.

if isfield(opts, 'casadi_path') && ~isempty(opts.casadi_path)
    local_progress_log(local_progress_enabled(opts), 'IPOPT', 'adding CasADi path: %s', opts.casadi_path);
    addpath(opts.casadi_path);
end
import casadi.*
local_progress_log(local_progress_enabled(opts), 'IPOPT', 'CasADi imported; nvars %d, constraints %d', numel(y0), numel(problem_definition.lbg));

% Build the shared IPOPT options used by both stages.
base_options = local_build_nlpsol_options(opts);
local_progress_log(local_progress_enabled(opts), 'IPOPT', 'options ready: tol %.3g, constr %.3g, print level %d', ...
    base_options.ipopt.tol, base_options.ipopt.constr_viol_tol, base_options.ipopt.print_level);

% Run the optional zero-objective feasibility stage.
run_stage1 = true;
if isfield(opts, 'run_stage1') && ~isempty(opts.run_stage1)
    run_stage1 = logical(opts.run_stage1);
end
if run_stage1
    stage1_nlp = problem_definition.nlp;
    stage1_nlp.f = 0 * problem_definition.nlp.f;
    stage1_options = base_options;
    stage1_options.ipopt.max_iter = opts.max_iter_stage1;
    local_progress_log(local_progress_enabled(opts), 'IPOPT', 'building stage 1 solver with max_iter %d', opts.max_iter_stage1);
    stage1_solver = nlpsol('rpo2ds_ipopt_stage1', 'ipopt', stage1_nlp, stage1_options);
    local_progress_log(local_progress_enabled(opts), 'IPOPT', 'starting stage 1 solve');
    stage1_result = stage1_solver('x0', y0(:), ...
        'lbx', problem_definition.lbx, ...
        'ubx', problem_definition.ubx, ...
        'lbg', problem_definition.lbg, ...
        'ubg', problem_definition.ubg);
    y_stage1 = full(stage1_result.x);
    output_stage1 = local_collect_output(stage1_solver, full(stage1_result.f), 'stage1');
    exitflag_stage1 = local_map_ipopt_status(output_stage1.return_status);
    local_progress_log(local_progress_enabled(opts), 'IPOPT', 'stage 1 complete: status %s, iter %s, objective %.6g', ...
        output_stage1.return_status, local_stats_value(output_stage1.stats, 'iter_count'), output_stage1.objective);
else
    local_progress_log(local_progress_enabled(opts), 'IPOPT', 'stage 1 skipped by opts.run_stage1');
    y_stage1 = y0(:);
    output_stage1 = struct('stage', 'stage1', ...
        'return_status', 'Stage 1 skipped by opts.run_stage1', ...
        'objective', 0, ...
        'stats', struct());
    exitflag_stage1 = NaN;
    stage1_options = base_options;
    stage1_options.ipopt.max_iter = opts.max_iter_stage1;
end

% Run the regularized delta-v objective stage from the selected warm start.
stage2_options = base_options;
stage2_options.ipopt.max_iter = opts.max_iter_stage2;
local_progress_log(local_progress_enabled(opts), 'IPOPT', 'building stage 2 solver with max_iter %d', opts.max_iter_stage2);
stage2_solver = nlpsol('rpo2ds_ipopt_stage2', 'ipopt', problem_definition.nlp, stage2_options);
local_progress_log(local_progress_enabled(opts), 'IPOPT', 'starting stage 2 solve');
stage2_result = stage2_solver('x0', y_stage1, ...
    'lbx', problem_definition.lbx, ...
    'ubx', problem_definition.ubx, ...
    'lbg', problem_definition.lbg, ...
    'ubg', problem_definition.ubg);
y_opt = full(stage2_result.x);
J_opt = full(stage2_result.f);
output_stage2 = local_collect_output(stage2_solver, J_opt, 'stage2');
exitflag_stage2 = local_map_ipopt_status(output_stage2.return_status);
local_progress_log(local_progress_enabled(opts), 'IPOPT', 'stage 2 complete: status %s, iter %s, objective %.6g', ...
    output_stage2.return_status, local_stats_value(output_stage2.stats, 'iter_count'), output_stage2.objective);

% Return the same solve payload shape expected by RPO2DS_transfer.
solve_out = struct();
solve_out.y_stage1 = y_stage1;
solve_out.y_opt = y_opt;
solve_out.J_opt = J_opt;
solve_out.exitflag = exitflag_stage2;
solve_out.output = output_stage2;
solve_out.output_stage1 = output_stage1;
solve_out.output_stage2 = output_stage2;
solve_out.exitflag_stage1 = exitflag_stage1;
solve_out.exitflag_stage2 = exitflag_stage2;
solve_out.options_stage1 = stage1_options;
solve_out.options_stage2 = stage2_options;
solve_out.options = stage2_options;
end

function log_config = local_progress_enabled(opts)
%LOCAL_PROGRESS_ENABLED Read the optional console progress-log flag from the opts struct.

log_config = false;
if isfield(opts, 'progress_logs') && ~isempty(opts.progress_logs)
    log_config = opts.progress_logs;
end
end

function value_text = local_stats_value(stats, field_name)
%LOCAL_STATS_VALUE Convert one optional IPOPT stats field into printable text.

if ~isstruct(stats) || ~isfield(stats, field_name)
    value_text = 'n/a';
    return
end
value = stats.(field_name);
if isnumeric(value) || islogical(value)
    value_text = num2str(value);
elseif isstring(value)
    value_text = char(value);
else
    value_text = char(string(value));
end
end

function local_progress_log(log_config, component, message, varargin)
%LOCAL_PROGRESS_LOG Print one timestamped progress message and flush the MATLAB console.

[enabled, t_start] = local_progress_state(log_config);
if ~enabled
    return
end
time_text = char(datetime('now', 'Format', 'HH:mm:ss'));
fprintf('[%s | +%.1fs] %s | %s\n', time_text, toc(t_start), component, sprintf(message, varargin{:}));
drawnow('limitrate');
end

function [enabled, t_start] = local_progress_state(log_config)
%LOCAL_PROGRESS_STATE Resolve boolean and elapsed-time origin for console progress logs.

enabled = false;
t_start = tic;
if isstruct(log_config)
    if isfield(log_config, 'enabled') && ~isempty(log_config.enabled)
        enabled = logical(log_config.enabled);
    end
    if isfield(log_config, 't_start') && ~isempty(log_config.t_start)
        t_start = log_config.t_start;
    end
else
    enabled = ~isempty(log_config) && logical(log_config);
end
if ~isstruct(log_config) && enabled && isappdata(0, 'RPO_progress_log_t_start')
    t_start = getappdata(0, 'RPO_progress_log_t_start');
end
end

function nlpsol_options = local_build_nlpsol_options(opts)
%LOCAL_BUILD_NLPSOL_OPTIONS Build the shared CasADi nlpsol and IPOPT options.

nlpsol_options = struct();
nlpsol_options.expand = false;
nlpsol_options.print_time = false;
nlpsol_options.show_eval_warnings = false;
nlpsol_options.ipopt = struct();
nlpsol_options.ipopt.tol = opts.optimality_tolerance;
nlpsol_options.ipopt.constr_viol_tol = opts.constraint_tolerance;
nlpsol_options.ipopt.acceptable_tol = max(opts.optimality_tolerance, 1e-6);
nlpsol_options.ipopt.acceptable_constr_viol_tol = max(opts.constraint_tolerance, 1e-8);
nlpsol_options.ipopt.acceptable_iter = 5;
nlpsol_options.ipopt.hessian_approximation = 'limited-memory';
nlpsol_options.ipopt.print_level = local_ipopt_print_level(opts);
end

function print_level = local_ipopt_print_level(opts)
%LOCAL_IPOPT_PRINT_LEVEL Convert the existing solver display string into an IPOPT print level.

print_level = 5;
if isfield(opts, 'solver_display') && strcmpi(string(opts.solver_display), "off")
    print_level = 0;
end
end

function output = local_collect_output(solver, objective, stage_name)
%LOCAL_COLLECT_OUTPUT Collect the raw IPOPT stats and compact fields used by reports.

stats = solver.stats();
return_status = '';
if isfield(stats, 'return_status')
    return_status = stats.return_status;
end
output = struct();
output.stage = stage_name;
output.return_status = return_status;
output.objective = objective;
output.stats = stats;
end

function exitflag = local_map_ipopt_status(return_status)
%LOCAL_MAP_IPOPT_STATUS Map IPOPT text statuses onto the numeric flags used by existing reports.

status = string(return_status);
if status == "Solve_Succeeded" || status == "Solved_To_Acceptable_Level"
    exitflag = 1;
elseif status == "Maximum_Iterations_Exceeded"
    exitflag = 0;
else
    exitflag = -1;
end
end
