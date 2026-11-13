function result = compute_manifolds(cfg,reference_direction)
%COMPUTE_MANIFOLDS Propagate both signs of the historical unstable velocity perturbation.
%
% PROTOTYPE:
%   result = sel2rpo.compute_manifolds(cfg, reference_direction)
%
% DESCRIPTION:
%   Recompute the reference halo and its monodromy matrix. Transport
%   the dominant unstable eigenvector along the orbit, normalize the
%   full six-component direction and perturb only its velocity part.
%   This reproduces the historical illustration, rather than displacing
%   the full state along a manifold eigenvector.
%
% INPUT:
%   cfg  Configuration with meta and sel2_to_rpo fields from
%        config/reproduction_config.m; physical units follow field suffixes.
%   cfg.sel2_to_rpo.manifolds sets phase count, perturbation and arc sampling.
%   reference_direction  [6 x 1] Frozen normalized-state direction used
%                          only to align the initial eigenvector sign.
%
% OUTPUT:
%   result  Result with data.halo_states [K+1 x 6], phase_times_nd
%             [K+1 x 1], periodic_orbit, constants and arcs.plus/minus
%             (K cells each). Each arc has xx_nd [M x 6], t_nd [M x 1]
%             and t_days. SE rotating states use [x,y,z,vx,vy,vz],
%             DU_SE/VU_SE and TU_SE; physical times use days. K and M
%             are phase_count and arc_samples. Metrics contain the
%             unstable multiplier [-] and normalized halo closure.

% Reconstruct the halo and integrate an identity-initialized STM over one period.
halo_result = sel2rpo.compute_halo(cfg);
halo = halo_result.data;
s = cfg.sel2_to_rpo.manifolds;
period = halo.solution.period_nd; initial = halo.solution.x0_nd;
options = odeset('RelTol',s.relative_tolerance,'AbsTol',s.absolute_tolerance);
identity = eye(6);
[time,extended] = ode113(@(t,x) rhs_SE3BP_3D_STM(t,x,halo.constants.mu_SE_nd), ...
    linspace(0,period,s.phase_count+1),[initial;identity(:)],options);

% Select the largest-magnitude multiplier and align its real direction.
monodromy = reshape(extended(end,7:42),6,6);
[vectors,values] = eig(monodromy);
[~,unstable] = max(abs(diag(values)));
direction = real(vectors(:,unstable)); direction = direction/norm(direction);
if dot(direction,reference_direction)<0, direction=-direction; end
arc_time = linspace(0,s.horizon_periods*period,s.arc_samples);
arcs.plus = cell(s.phase_count,1); arcs.minus = cell(s.phase_count,1);

% Transport the direction and propagate both velocity-only perturbations.
for phase = 1:s.phase_count
    transition = reshape(extended(phase,7:42),6,6);
    transported = transition*direction; transported = transported/norm(transported);
    for sign = [1,-1]
        state = extended(phase,1:6).';
        state(4:6) = state(4:6)+sign*s.velocity_perturbation*transported(4:6);
        [t,x] = ode113(@(t,x) rhs_SE3BP_3D(t,x,halo.constants.mu_SE_nd),arc_time,state,options);
        arc = struct('t_nd',t,'t_days',t*halo.constants.TU_SE_days,'xx_nd',x);
        if sign==1, arcs.plus{phase}=arc; else, arcs.minus{phase}=arc; end
    end
    if mod(phase,10)==0, fprintf('Unstable branches: %d/%d phases\n',phase,s.phase_count); end
end

% Store both branch illustrations and the full-period closure diagnostic.
data = struct('constants',halo.constants, ...
    'periodic_orbit',struct('x0_nd',initial,'period_nd',period), ...
    'halo_states',extended(:,1:6),'phase_times_nd',time, ...
    'arcs',arcs,'unstable_multiplier',values(unstable,unstable), ...
    'reference_direction',direction);
metrics = struct('unstable_multiplier',data.unstable_multiplier, ...
    'closure_inf_nd',norm(extended(end,1:6).'-initial,inf));
result = rpo.finalize_result(struct('data',data,'metrics',metrics),cfg, ...
    'run_04_sel2_to_rpo: unstable directions','analytically reconstructed halo','generated_result');
end
