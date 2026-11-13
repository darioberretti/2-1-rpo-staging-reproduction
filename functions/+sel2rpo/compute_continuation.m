function result = compute_continuation(cfg,inputs,result_file)
%COMPUTE_CONTINUATION Reproduce the historical Sun-phase campaigns and merge their cases.
%
% PROTOTYPE:
%   result = sel2rpo.compute_continuation(cfg, inputs, result_file)
%
% DESCRIPTION:
%   Solve the zero/pi seed families, the later 1.3-rad family and
%   its positive extension, preserving historical phase grids and warm
%   starts. Each branch advances serially even when solver flags are
%   unfavorable. After each family, optionally save an intermediate
%   merge with status.completed=false; only the final merge is complete.
%   These saves record progress but do not provide checkpoint/resume.
%
% INPUT:
%   cfg  Configuration with meta and sel2_to_rpo fields from
%        config/reproduction_config.m; physical units follow field suffixes.
%   cfg.sel2_to_rpo.optimization and continuation set solve/grid settings.
%   inputs       Frozen inputs with three problems (y0, nlp, constants,
%                  label and meta.source_thetaSun0_rad [rad]) from
%                  data/inputs/run_04_reference_inputs.mat.
%   result_file  Destination MAT path, or empty text to disable saves.
%
% OUTPUT:
%   result       Provenance-bearing result with data.cases (merged
%                  solve_case records), problems and four campaign arrays.
%                  Metrics include case_count, runtime_seconds [s],
%                  delta_v_range_m_s [m/s], flight/elapsed_range_days
%                  [days], and maximum normalized equality residual.
%                  Phase duplicates use circular tolerance; lower cost
%                  wins, with later campaigns winning equal-cost ties.

s = cfg.sel2_to_rpo.optimization; grid = cfg.sel2_to_rpo.continuation;
campaigns = cell(4,1); clock_start = tic;
fprintf('Continuation starts from frozen unoptimized nodes for all three seed problems.\n');

% Solve each seed family with separate negative and positive warm starts.
for family=1:3
    problem = inputs.problems(family);
    theta = problem.meta.source_thetaSun0_rad;
    center = sel2rpo.solve_case(problem,problem.y0,theta,s,s.constraint_tolerance, ...
        family,problem.label,0,1);
    % Each historical branch starts from the original seed, not the separately
    % re-solved grid center. They differ for the later 1.3-rad alignment.
    branch_initial_y = center.y_sol;

    % Align only the later family center to the historical phase grid.
    if family==3
        reference_grid = (0:255)*grid.phase_step_rad;
        distance = abs(atan2(sin(reference_grid-theta),cos(reference_grid-theta)));
        [~,nearest] = min(distance); theta = reference_grid(nearest);
        center = sel2rpo.solve_case(problem,center.y_sol,theta,s,s.constraint_tolerance, ...
            family,problem.label,0,1);
    end

    % Solve both branches from the retained initial vector, then save progress.
    negative = run_branch(problem,center,branch_initial_y,-1,grid,s,s.constraint_tolerance,family,problem.label);
    positive = run_branch(problem,center,branch_initial_y,1,grid,s,s.constraint_tolerance,family,problem.label);
    campaigns{family} = [center;negative;positive];
    if family==3
        % Historical extraction selected positive case 128, before its tiny
        % floating-point duplicate at the clamped endpoint.
        extension_seed = positive(128);
    end
    result = package(cfg,inputs,campaigns(1:family),grid.merge_tolerance_rad,toc(clock_start),false);
    if ~isempty(result_file), rpo.save_result(result,result_file); end
end

% Extend the selected positive case with the historical relaxed tolerance.
problem = inputs.problems(3);
tolerance = s.extension_constraint_tolerance;
center = sel2rpo.solve_case(problem,extension_seed.y_sol,extension_seed.theta_sun0_rad, ...
    s,tolerance,3,'positive_extension',0,1);
positive = run_branch(problem,center,extension_seed.y_sol,1,grid,s,tolerance,3,'positive_extension');
campaigns{4} = [center;positive];
result = package(cfg,inputs,campaigns,grid.merge_tolerance_rad,toc(clock_start),true);
if ~isempty(result_file), rpo.save_result(result,result_file); end
end

function cases = run_branch(problem,center,initial_y,direction,grid,s,tolerance,problem_index,label)
%RUN_BRANCH Advance a serial warm-start branch on the historical clamped phase grid.
%
% PROTOTYPE:
%   cases = run_branch(problem, center, initial_y, direction, grid, s, ...
%       tolerance, problem_index, label)
%
% DESCRIPTION:
%   Use repeated addition, clamp the final phase to the half-span
%   endpoint and pass every returned y_sol to the next case, regardless
%   of solver flags. The center is excluded from the returned array.
%
% INPUT:
%   problem, s   Frozen problem and optimization settings for solve_case.
%   center       solve_case record setting the initial Sun phase [rad].
%   initial_y    Initial branch decision vector, independent of center.y_sol.
%   direction    -1 or +1 for descending/ascending phases.
%   grid         phase_step_rad and half_span_rad [rad].
%   tolerance    Constraint tolerance used for this branch.
%   problem_index, label  Seed-family index and text label.
%
% OUTPUT:
%   cases        Column struct array of solve_case records in solve order.

theta = center.theta_sun0_rad;
limit = theta+direction*grid.half_span_rad;
y = initial_y; index=0;
while true
    next = theta+direction*grid.phase_step_rad;
    done = direction*(next-limit)>=0;
    if done, next=limit; end
    index=index+1;
    solved = sel2rpo.solve_case(problem,y,next,s,tolerance,problem_index,label,direction,index);
    if index==1, cases=solved; else, cases(end+1,1)=solved; end %#ok<AGROW>
    % Advance even when flags are unfavorable, matching the archived method.
    y=solved.y_sol; theta=next;
    if done, break; end
end
end

function result = package(cfg,inputs,campaigns,tolerance,runtime,completed)
%PACKAGE Merge phase duplicates and attach campaign provenance and metrics.
%
% PROTOTYPE:
%   result = package(cfg, inputs, campaigns, tolerance, runtime, completed)
%
% DESCRIPTION:
%   Visit later campaigns first, wrap phases to [0,2*pi) and keep
%   the lower-cost case within circular tolerance. Equal costs retain
%   the first visited case; finally sort by wrapped Sun phase.
%
% INPUT:
%   cfg, inputs  Scientific configuration and frozen seed problems.
%   campaigns    Cell array of solve_case struct arrays.
%   tolerance    Circular phase-merge tolerance [rad].
%   runtime      Cumulative elapsed runtime [s].
%   completed    Logical marking a final rather than intermediate merge.
%
% OUTPUT:
%   result       Merged cases, problems, original campaigns, cost/time
%                  ranges, residual/runtime metrics and status.completed.

all_cases = vertcat(campaigns{end:-1:1});
merged = all_cases([]);
for index=1:numel(all_cases)
    candidate=all_cases(index); candidate.theta_sun0_rad=mod(candidate.theta_sun0_rad,2*pi);
    phases=[merged.theta_sun0_rad];
    distance=abs(atan2(sin(phases-candidate.theta_sun0_rad),cos(phases-candidate.theta_sun0_rad)));
    match=find(distance<=tolerance,1);
    if isempty(match), merged(end+1,1)=candidate; %#ok<AGROW>
    elseif candidate.delta_v_m_s<merged(match).delta_v_m_s, merged(match)=candidate; end
end

% Sort the merged phases and summarize every retained case, including its residual.
[~,order]=sort([merged.theta_sun0_rad]); merged=merged(order);
data = struct('cases',merged,'problems',inputs.problems,'campaigns',{campaigns});
metrics = struct('case_count',numel(merged),'runtime_seconds',runtime, ...
    'delta_v_range_m_s',[min([merged.delta_v_m_s]),max([merged.delta_v_m_s])], ...
    'flight_range_days',[min([merged.flight_days]),max([merged.flight_days])], ...
    'elapsed_range_days',[min([merged.elapsed_days]),max([merged.elapsed_days])], ...
    'maximum_equality_residual',max([merged.equality_residual]));
result=rpo.finalize_result(struct('data',data,'metrics',metrics),cfg, ...
    'run_04_sel2_to_rpo: continuation','data/inputs/run_04_reference_inputs.mat: unoptimized seed problems','generated_result');
result.status.completed=completed;
end
