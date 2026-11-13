function trajectory = reconstruct_transfer(y, problem)
%RECONSTRUCT_TRANSFER Propagate shooting arcs independently in EM and SE rotating frames.
%
% PROTOTYPE:
%   trajectory = sel2rpo.reconstruct_transfer(y, problem)
%
% DESCRIPTION:
%   Decode the historical case-2a regularized impulses and boundary
%   phases. Interior states are post-impulse. Reconstruct boundaries
%   using the historical EMB-centered SE secondary, propagate each
%   BCR4BP arc, and measure its endpoint defect without joining arcs.
%
% INPUT:
%   y        [10*N-11 x 1] Decision vector for N shooting nodes:
%            U1, interior (Xj,Uj), UN, N-1 segment durations, SEL2
%            phase and RPO phase. Xj = [x,y,z,vx,vy,vz] uses the
%            EM rotating frame, DU_EM and VU_EM. Uj = [u,w,s]
%            defines a regularized impulse; durations/RPO phase
%            use TU_EM and SEL2 phase uses TU_SE.
%   problem  Frozen shooting problem from inputs.problems in
%            data/inputs/run_04_reference_inputs.mat. Contains nlp
%            (node maps, boundary seeds, phases, model constants,
%            constraints and ODE tolerances) and constants
%            (TU_days [days/TU_EM], VU_km_s [km/s per VU_EM]).
%
% OUTPUT:
%   trajectory  Struct with N node rows in nodes_em/nodes_se [N x 6],
%                 node_times_em [N x 1, TU_EM], impulses_em [N x 3, VU_EM]
%                 and N-1 segment cells containing time_em [M x 1],
%                 states_em/states_se [M x 6]. M varies with ODE output.
%                 SE states use DU_SE/VU_SE; EM states use DU_EM/VU_EM.
%                 defects [N-1 x 6] stores arrival minus next-node
%                 position, and next-node minus arrival minus impulse
%                 velocity, in EM units. Also returns halo_states_se,
%                 waiting/flight/elapsed_days [days], total_delta_v_m_s
%                 [m/s], halo_phase [TU_SE], rpo_phase [TU_EM],
%                 theta_sun0_rad [rad] and both mass ratios [-].

% Decode node impulses and interior post-impulse states.
nlp = problem.nlp;
node_count = nlp.N;
impulses = zeros(node_count,3);
nodes = zeros(node_count,6);
for node = 1:node_count
    regularized = y(nlp.map.node(node).U);
    u = regularized(1); w = regularized(2); s = regularized(3);
    impulses(node,:) = [u^2-w^2-s^2, 2*u*w, 2*u*s];
    state_indices = nlp.map.node(node).X;
    if ~isempty(state_indices), nodes(node,:) = y(state_indices).'; end
end

% Convert the wrapped halo wait to the EM timeline before adding flight times.
halo_phase = y(nlp.map.theta_sel2);
waiting_se = mod(halo_phase-nlp.theta_halo0_nd, nlp.sel2_period_nd);
waiting_em = waiting_se*nlp.omega_em_rad_s/nlp.omega_se_rad_s;
durations = y(nlp.map.dt_nodes);
node_times = waiting_em + [0;cumsum(durations(:))];

% Reconstruct the SE halo departure and CR3BP RPO arrival boundaries.
boundary_options = odeset('RelTol',1e-12,'AbsTol',1e-14);
if halo_phase == nlp.theta_halo0_nd
    halo_state = nlp.sel2_traj_phase.x0_se_nd(:);
else
    halo_solution = ode113(@(t,x) rhs_SE3BP_3D(t,x,nlp.mu_se_nd), ...
        [nlp.theta_halo0_nd,halo_phase], nlp.sel2_traj_phase.x0_se_nd(:), boundary_options);
    halo_state = halo_solution.y(:,end);
end
sun_phase = nlp.theta_sun0_rad + nlp.ws_sun_em_nd*node_times(1);
nodes(1,:) = se_cr3bp_to_emb_bcr4bp_state(halo_state.',sun_phase, ...
    'alpha_plus_theta_sun_const_rad',nlp.alpha_plus_theta_const, ...
    'l_se_km',nlp.l_se_km,'l_em_km',nlp.l_em_km, ...
    'omega_se_rad_s',nlp.omega_se_rad_s,'omega_em_rad_s',nlp.omega_em_rad_s, ...
    'mu_se_nd',nlp.mu_se_nd,'se_secondary_center','emb');
nodes(1,4:6) = nodes(1,4:6) + impulses(1,:);
rpo_seed = nlp.rpo_traj_phase.x0_rpo_planar(:);
rpo_state = [rpo_seed(1:2);0;rpo_seed(3:4);0];
nodes(end,:) = propagate(rpo_state,0,y(nlp.map.theta_rpo), ...
    nlp.mu_em,'3D','CR3BP','false').';

% Propagate each EM shooting arc separately and retain its continuity defect.
segments = cell(node_count-1,1);
defects = zeros(node_count-1,6);
for node = 1:node_count-1
    [arrival,~,~,states,times] = propagate(nodes(node,:).', ...
        node_times(node),node_times(node+1),nlp.mu_em,'3D','BCR4BP','false', ...
        [],[],nlp.theta_sun0_rad,nlp.abs_tol,nlp.rel_tol);
    sun_phases = nlp.theta_sun0_rad + nlp.ws_sun_em_nd*times;
    states_se = convert_to_se(states(:,1:6),sun_phases,nlp);
    segments{node} = struct('time_em',times,'states_em',states(:,1:6),'states_se',states_se);
    defects(node,1:3) = arrival(1:3).' - nodes(node+1,1:3);
    defects(node,4:6) = nodes(node+1,4:6) - arrival(4:6).' - impulses(node+1,:);
end

% Convert nodes and sample the reference halo for the two-frame plots.
node_sun_phases = nlp.theta_sun0_rad + nlp.ws_sun_em_nd*node_times;
nodes_se = convert_to_se(nodes,node_sun_phases,nlp);
halo_grid = linspace(0,nlp.sel2_period_nd,500);
[~,halo_states] = ode113(@(t,x) rhs_SE3BP_3D(t,x,nlp.mu_se_nd), ...
    halo_grid,nlp.sel2_traj_phase.x0_se_nd(:),boundary_options);

% Return normalized geometry together with physical cost and time metrics.
trajectory = struct('segments',{segments},'nodes_em',nodes,'nodes_se',nodes_se, ...
    'node_times_em',node_times,'impulses_em',impulses,'defects',defects, ...
    'halo_states_se',halo_states,'waiting_days',waiting_em*problem.constants.TU_days, ...
    'flight_days',sum(durations)*problem.constants.TU_days, ...
    'elapsed_days',node_times(end)*problem.constants.TU_days, ...
    'total_delta_v_m_s',sum(vecnorm(impulses,2,2))*problem.constants.VU_km_s*1000, ...
    'halo_phase',halo_phase,'rpo_phase',y(nlp.map.theta_rpo), ...
    'theta_sun0_rad',nlp.theta_sun0_rad,'mu_em',nlp.mu_em,'mu_se',nlp.mu_se_nd);
end

function states_se = convert_to_se(states_em,sun_phases,nlp)
%CONVERT_TO_SE Convert EM state rows to SE units with the historical EMB secondary.
% states_em/states_se are [M x 6]; sun_phases are scalar or M angles [rad].

states_se = emb_bcr4bp_to_se_cr3bp_state(states_em,sun_phases, ...
    'alpha_plus_theta_sun_const_rad',nlp.alpha_plus_theta_const, ...
    'l_se_km',nlp.l_se_km,'l_em_km',nlp.l_em_km, ...
    'omega_se_rad_s',nlp.omega_se_rad_s,'omega_em_rad_s',nlp.omega_em_rad_s, ...
    'mu_se_nd',nlp.mu_se_nd,'se_secondary_center','emb');
end
