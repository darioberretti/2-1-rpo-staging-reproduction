function result = compute_halo(cfg)
%COMPUTE_HALO Correct and continue a Richardson seed to the target SEL2 halo.
%
% PROTOTYPE:
%   result = sel2rpo.compute_halo(cfg)
%
% DESCRIPTION:
%   Start from a third-order L2 seed, correct x0, vy0 and half-period
%   at fixed z0, then scale z0 until the target amplitude is bracketed.
%   Bisect z0 with the historical warm starts. Independently propagate
%   the final orbit to check closure and sampled vertical amplitude.
%
% INPUT:
%   cfg  Configuration with meta and sel2_to_rpo fields from
%        config/reproduction_config.m; physical units follow field suffixes.
%   cfg.sel2_to_rpo.halo contains constants, targets and solver settings.
%
% OUTPUT:
%   result  Provenance-bearing result with data.constants, targets,
%             solver, seed, continuation_history, solution and trajectory.
%             trajectory.xx_nd [M x 6] is [x,y,z,vx,vy,vz] in the SE
%             rotating barycentric frame (DU_SE, VU_SE); t_nd [M x 1]
%             uses TU_SE and t_days uses days. solution/metrics report
%             period in TU_SE/days, Ay/Az in km and closure in normalized
%             state components. M is the configured final sample count.

settings = cfg.sel2_to_rpo.halo;
c = settings.constants; target = settings.targets; solver = settings.solver;

% Build the analytical seed in the historical SE normalization.
seed = richardson_seed_halo_se3bp(c.mu_SE_nd,'L2',target.Az_target_nd, ...
    'north_south',target.north_south_sign,'t_nd',0);
state = seed.x0_nd(:); half_period = 0.5*seed.T_est_nd;
history = struct([]);

% Continue corrected seeds in z0, then bisect the amplitude bracket.
for step = 1:solver.continuation_max_steps
    corrected = correct(state,half_period,c,solver);
    if step==1, history=corrected; else, history(end+1)=corrected; end %#ok<AGROW>
    fprintf('Halo continuation %d: Az = %.6f km\n',step,corrected.Az_km);
    if abs(corrected.Az_km-target.Az_target_km)<=solver.Az_tolerance_km
        break;
    elseif corrected.Az_km<target.Az_target_km
        % Reuse the above-target corrected state, as in the original recipe.
        above = history(end-1); below = corrected;
        high_z = above.seed(3); low_z = below.seed(3);
        state = above.state; half_period = above.period/2;
        for iteration = 1:solver.bisection_max_steps
            mid_z = (high_z+low_z)/2; state(3) = mid_z;
            corrected = correct(state,half_period,c,solver);
            history(end+1) = corrected; %#ok<AGROW>
            fprintf('Halo bisection %d: Az = %.6f km\n',iteration,corrected.Az_km);
            if abs(corrected.Az_km-target.Az_target_km)<=solver.Az_tolerance_km
                break;
            end
            if corrected.Az_km<target.Az_target_km, low_z=mid_z; else, high_z=mid_z; end
            state = corrected.state; half_period = max(corrected.period/2,solver.dc_min_half_period_nd);
        end
        break;
    end
    state = corrected.state;
    state(3) = solver.continuation_z0_scale*state(3);
    half_period = max(corrected.period/2,solver.dc_min_half_period_nd);
end

% Verify the final orbit independently at the reporting resolution.
options = odeset('RelTol',solver.ode_rel_tol,'AbsTol',solver.ode_abs_tol);
[time,states] = ode113(@(t,x) rhs_SE3BP_3D(t,x,c.mu_SE_nd), ...
    linspace(0,corrected.period,solver.n_output_samples_final),corrected.state,options);
closure = norm(states(end,:).'-corrected.state,inf);
az = max(abs(states(:,3)))*c.L_SE_km;
if closure>1e-8 || abs(az-target.Az_target_km)>solver.Az_tolerance_km
    error('sel2rpo:HaloCorrection','Halo correction did not meet closure/amplitude tolerances.');
end

% Package physical amplitudes, normalized closure and reconstruction history.
solution = struct('x0_nd',corrected.state,'period_nd',corrected.period, ...
    'period_days',corrected.period*c.TU_SE_days, ...
    'Ay_achieved_km',max(abs(states(:,2)))*c.L_SE_km,'Az_achieved_km',az, ...
    'closure_inf_nd',closure,'continuation_steps_done',numel(history));
data = struct('constants',c,'targets',target,'solver',solver,'seed',seed, ...
    'continuation_history',history,'solution',solution, ...
    'trajectory',struct('t_nd',time,'t_days',time*c.TU_SE_days,'xx_nd',states));
result = rpo.finalize_result(struct('data',data,'metrics',solution),cfg, ...
    'run_04_sel2_to_rpo: halo','analytical Richardson seed','generated_result');
end

function step = correct(state,half_period,c,solver)
%CORRECT Correct one fixed-z seed and measure its full-orbit amplitude.
%
% PROTOTYPE:
%   step = correct(state, half_period, c, solver)
%
% DESCRIPTION:
%   Reject a failed half-period correction, then sample one full
%   period with the configured continuation resolution.
%
% INPUT:
%   state        [6 x 1] SE rotating seed [x;y;z;vx;vy;vz], normalized.
%   half_period  Positive initial half-period guess [TU_SE].
%   c            Halo constants, including mu_SE_nd [-] and L_SE_km.
%   solver       cfg.sel2_to_rpo.halo.solver settings.
%
% OUTPUT:
%   step         Struct with seed, corrected state [6 x 1], period
%                  [TU_SE], sampled Az_km [km] and correction diagnostics.

dc = diffcorr_halo_halfperiod_se3bp(state,c.mu_SE_nd, ...
    'half_period_guess_nd',half_period,'max_iter',solver.dc_max_iter, ...
    'tol_res_inf',solver.dc_tol_res_inf_nd,'tol_step_inf',solver.dc_tol_step_inf_nd, ...
    'damping',solver.dc_damping,'use_pinv',solver.dc_use_pinv, ...
    'cond_rcond_min',solver.dc_cond_rcond_min,'ode_rel_tol',solver.ode_rel_tol, ...
    'ode_abs_tol',solver.ode_abs_tol,'min_half_period_nd',solver.dc_min_half_period_nd, ...
    'verbose',solver.dc_verbose);
if ~dc.success, error('sel2rpo:HaloCorrection','%s',dc.message); end
[~,states] = ode113(@(t,x) rhs_SE3BP_3D(t,x,c.mu_SE_nd), ...
    linspace(0,dc.T_nd,solver.n_output_samples_continuation),dc.x0_corrected_nd, ...
    odeset('RelTol',solver.ode_rel_tol,'AbsTol',solver.ode_abs_tol));
step = struct('seed',state,'state',dc.x0_corrected_nd,'period',dc.T_nd, ...
    'Az_km',max(abs(states(:,3)))*c.L_SE_km,'correction',dc);
end
