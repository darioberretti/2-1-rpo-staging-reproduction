function [vInfMag, dv_opt, tof_opt, y_opt, J_opt, exitflag, output, h] = RPO2DS_reachability(vInfPlusEarth, dv_budget, RPOMoonDatabase, xx0_RPO, params, graph, gradient_check, initial_guess)
%RPO2DS_REACHABILITY Run the historical two-stage fixed-budget reachability solve.
%
% PROTOTYPE:
%   [vInfMag, dv_opt, tof_opt, y_opt, J_opt, exitflag, output, h] = ...
%       RPO2DS_reachability(vInfPlusEarth, dv_budget, RPOMoonDatabase, ...
%       xx0_RPO, params, graph, gradient_check, initial_guess)
%
% DESCRIPTION:
%   Reconstruct the patched trajectory, place shooting nodes, regularize
%   the impulses, and run the original fmincon feasibility/speed objectives.
%   Run02 supplies archived patched starts through initial_guess; this skips
%   preprocessing but still executes the optimizations. The backend has its
%   own propagator/constraints selected temporarily by compute_reachability.
%   See docs/REPRODUCIBILITY.md (Run 02) for fresh-run agreement limits.
%
% INPUT:
%   vInfPlusEarth   [3 elements] Earth-relative target direction/seed [km/s].
%   dv_budget       Total normalized impulse budget [VU].
%   RPOMoonDatabase  Historical dense or phase-aware RPO-Moon table.
%   xx0_RPO         [4 elements] [x, y, vx, vy] in rotating DU/VU.
%   params          [node_count, minimum_Earth_altitude_km].
%   graph           Optional plotting flag, default 0.
%   gradient_check  Optional derivative-check flag, default 0.
%   initial_guess   Optional struct with full_solution patch table and
%                     reported_altitude_km; absent/empty rebuilds the patch.
%
% OUTPUT:
%   vInfMag         Optimized outgoing Earth-relative speed [km/s].
%   dv_opt          Total regularized impulse magnitude [m/s].
%   tof_opt         Optimized transfer flight time [days].
%   y_opt           Historical [10*N + 1 x 1] decision vector: node blocks
%                     [X(6); U(3)], N absolute times [TU], and vInfMag [km/s].
%   J_opt           Final negative-speed objective; exitflag/output are
%                     the final fmincon status and diagnostics.
%   h               Earth-perigee altitude reported by the patch search [km].

if nargin < 6
    graph = 0;
end

if nargin < 7
    gradient_check = 0;
end
n_nodes = params(1);
hStart = params(2);
% Constants
mu_EM = 0.01215;                     % Earth–Moon system mass parameter [nondimensional]
muE = 398600;                        % Earth gravitational parameter [km^3/s^2]
DU = 384400;                         % Distance between Earth and Moon [km]
om_EM = 2.665e-6;                    % Angular velocity of Earth and Moon around EM barycenter [rad/s]
TU = 1 / om_EM;                      % Time unit in Earth–Moon system [s]
VU = DU / TU;                        % Velocity unit in Earth–Moon system [km/s]
Re = 6378;                           % Earth radius [km]
Rm = 1738;                           % Moon radius [km]
rSOI = (mu_EM / (1 - mu_EM))^(2/5); % Moon's sphere of influence (SOI) radius in EM system (Vallado) [nondimensional]
hminEarth = hStart;                     % Minimum altitude from Earth's surface [km]
hminMoon = 100;                      % Minimum altitude from Moon's surface [km]
rminEarth = Re + hminEarth;          % Minimum distance from Earth's center [km]
rminMoon = Rm + hminMoon;            % Minimum distance from Moon's center [km]
rminEarth_DU = rminEarth / DU;       % Minimum distance from Earth's center in adimensional units[DU]
rminMoon_DU = rminMoon / DU;         % Minimum distance from Moon's center in adimensional units[DU]
maxTOF = 100;                        % Days
maxTOF = maxTOF * 86400 / TU;        % TU

% Build or restore the historical pre-optimization patch
if nargin < 8 || isempty(initial_guess)
    patch_parameters = struct('N', 101, 'hStart', hStart, ...
        'sweep_rad', pi, 'nSol', 1, 'feasibleLimit', 1e-2, ...
        'nMoonOrbitSamples', 100, 'time_mode', 'optimization_only', ...
        'maxTOF', 100, 'geometry_mode', 'free_geometry', ...
        'perigee_sampling_type', 1);
    [FullSolution, ~, h] = RPO_Moon_Earth_guess( ...
        vInfPlusEarth, RPOMoonDatabase, patch_parameters);
else
    FullSolution = initial_guess.full_solution;
    h = initial_guess.reported_altitude_km;
end

% Make the trajectory ready for optimization
% Flyby insertion
v_in = FullSolution.vInfMinusMoon;
v_out = FullSolution.vInfPlusMoon;
posMoon = FullSolution.RI;
[Rp, Vp] = get_perilune_state(v_in, v_out);
RpECI = RTH2ECI_moon(posMoon, Rp, [], 1);
VpECI = RTH2ECI_moon(posMoon, Vp, [], 0);
% Propagate forward and backward the perilune state until it crosses the
% Moon SOI (in the RTH frame)
x0 = [RpECI(:)/DU; VpECI(:)/VU];
time = FullSolution.trueAnomalyMoonFlyby;
xx0_peri = in2rot(x0',time,mu_EM);
tspanF  = [time, time+1];
tspanB  = [time, time-1];
opts   = odeset('RelTol',1e-10,'AbsTol',1e-12, ...
                'Events', @(t,x) soi_event_CR3BP(t,x,rSOI,mu_EM));
[tF, xxF] = ode113(@(t,x) rhs_CR3BP_3D(t,x,mu_EM), tspanF, xx0_peri, opts);
[tB, xxB] = ode113(@(t,x) rhs_CR3BP_3D(t,x,mu_EM), tspanB, xx0_peri, opts);
xxFB = [xxB(end:-1:2,:); xxF];
ttFB = [tB(end:-1:2,:); tF];
% Total flyby time
tIncFlyby = tB(1) - tB(end);
tTotFlyby = ttFB(end) - ttFB(1);

% CR3BP propagation (RPO-Moon)
xx0RPO = [xx0_RPO(1), xx0_RPO(2), 0, xx0_RPO(3), xx0_RPO(4), 0]';
if ismember("orbitPhase", RPOMoonDatabase.Properties.VariableNames)
    phase = FullSolution.phase;
    [xx0RPO]  = propagate(xx0RPO,0,phase,mu_EM,'3D','CR3BP','false');
end
x0 = xx0RPO;
x0(4:6) = x0(4:6) + FullSolution.dv0CR3BP';
% Adjusting time of departure to match true anomaly of the Moon on the
% Moon-Earth leg
angle = ttFB(1);
angleDifference = angdiff(FullSolution.trueAnomalyMoonRPO, angle);
t0 = angleDifference;
tf = FullSolution.T_RPOMoon + t0;
[~,~,~, xx, tt]  = propagate(x0,t0,tf,mu_EM,'3D','CR3BP','false');

% True time flyby
ttFB_from_zero = ttFB - ttFB(1);
ttFB_true = ttFB_from_zero + tf;

% TBP propagation (Moon-Earth leg)
rhs = @(t,y) ode_2bp(t,y,muE);
xx0 = [FullSolution.RI'; FullSolution.VI'];
t0_lamb = (tt(end) + tIncFlyby) * TU;
tf_lamb = FullSolution.T_lambert + t0_lamb;
[t, y] = ode113(rhs, [t0_lamb tf_lamb], xx0, opts);
% Keep only after SOI
tSOI = (tt(end) + tTotFlyby);
t = t/TU;
indices = find(t > tSOI);
y = y(indices, :);
t = t(indices);
y_adim(:,1:3) = y(:,1:3) / DU;
y_adim(:,4:6) = y(:,4:6) / VU;
xxRot = in2rot(y_adim, t, mu_EM);

% Put together times and vectors
tTot = [tt; ttFB_true(2:end); t]; 
xxTot = [xx; xxFB(2:end,:); xxRot];
% Add DVp to the last state (for fmincon)
vp = y(end,4:6);
rp = y(end,1:3);
vpaft = vp+FullSolution.deltavp;
xxpaft_rot = in2rot([rp/DU, vpaft/VU], t(end),mu_EM);
xxTot(end,:) = xxpaft_rot;
% Pass to the next phase the index of the RPO-Moon phase
RPO_idx = numel(tt);
moon_peri_idx = RPO_idx + numel(tB) - 2; % -2 for lost indexes in merging

% Check for duplicate values in tTot
[unique_tTot, idx] = unique(tTot);
if numel(tTot) ~= numel(unique_tTot)
    fprintf('Warning: tTot contains duplicate values.\n');
    duplicate_indices = setdiff(1:numel(tTot), unique(idx));
    fprintf('Indices of duplicate values: %s\n', mat2str(duplicate_indices));
else
    fprintf('tTot does not contain duplicate values.\n');
end
if graph
% Moon swing by CR3BP
figure;
hold on;
grid on;
plot3(xxF(1,1), xxF(1,2), xxF(1,3), 'r.', 'MarkerSize', 20, 'DisplayName','Perilune');
plot3(1-mu_EM, 0, 0, 'k.', 'MarkerSize', 20, 'DisplayName','Moon');
plot3(xxF(:,1), xxF(:,2), xxF(:,3), 'r', 'LineWidth', 2, 'DisplayName','Forward');
plot3(xxB(:,1), xxB(:,2), xxB(:,3), 'b', 'LineWidth', 2, 'DisplayName','Backward');
xlabel('X (DU)', 'Interpreter', 'latex');
ylabel('Y (DU)', 'Interpreter', 'latex');
zlabel('Z (DU)', 'Interpreter', 'latex');
title('Moon swing-by - Adimensional EMB Centered Rotating Frame', 'Interpreter','latex');
set(gcf,"Theme","light");
legend('show', 'interpreter', 'latex');
axis equal;
xlim([0.8 1.05]);
ylim([min(xxFB(:,2)) max(xxFB(:,2))]);
grid on;
hold off;
set(gca,'FontSize',16,'FontName','times new roman');

% Full trajectory CR3BP
figure;
hold on;
plot3(xx(:,1),xx(:,2),xx(:,3), 'LineWidth',1.5,  'DisplayName', 'RPO-Moon leg');
plot3(xxRot(:,1), xxRot(:,2),xxRot(:,3), 'LineWidth',1.5,  'DisplayName', 'Moon-Earth leg');
theta = linspace(0, 2*pi, 400);
xc = (1 - mu_EM) + rSOI*cos(theta);
yc = rSOI*sin(theta);
zc = - ones(size(theta));
fill3(xc, yc, zc, [1.0  1.0  0.6], 'FaceAlpha', 1, 'EdgeColor', 'none', 'DisplayName', 'Moon SOI');
plot3(xxFB(:,1), xxFB(:,2), xxFB(:,3), 'LineWidth', 1.5, 'DisplayName','Moon Flyby leg');
plot3(xx(1,1),xx(1,2),xx(1,3), 'o', 'Color','b','MarkerSize', 10, 'LineWidth', 3, 'DisplayName','Initial point');
plot3(xxRot(end,1), xxRot(end,2),xxRot(end,3), 'o', 'Color','r','MarkerSize', 10, 'LineWidth', 3, 'DisplayName','Final point');
axis equal;
xlabel('X (DU)', 'Interpreter', 'latex');
ylabel('Y (DU)', 'Interpreter', 'latex');
zlabel('Z (DU)', 'Interpreter', 'latex');
title('Initial Guess Trajectory - Adimensional EMB Centered Rotating Frame', 'Interpreter','latex')
set(gcf,"Theme","light");
legend('show', 'Interpreter', 'Latex');
xlim([min(xxTot(:,1))*1.05 max(xxTot(:,1))*1.75]);
ylim([min(xxTot(:,2)) max(xxTot(:,2))]);
grid on;
hold off;
set(gca,'FontSize',16,'FontName','times new roman');
end
% Optimization
% Node generation and initial guess
% n_nodes = 40;  % Total number of nodes

% Create node times
density_sep_idx = round(RPO_idx * 0.99); % after this index, the node density increases
dense_nodes_multiplier = 0.15; % multiplier for total nodes assigned to the dense part
dense_nodes = round(n_nodes * dense_nodes_multiplier);
nondense_nodes = n_nodes - dense_nodes;
t_nodes_IG_nondense = linspace(tTot(1), tTot(density_sep_idx), nondense_nodes).';  % column vector
t_nodes_IG_dense = linspace(tTot(density_sep_idx), tTot(end), dense_nodes + 1).';  % column vector
t_nodes_IG = [t_nodes_IG_nondense; t_nodes_IG_dense(2:end)];

% Interpolate states at node times using pchip
X_nodes_IG = interp1(tTot, xxTot, t_nodes_IG, 'pchip');

% Add perilune node
t_peri = tTot(moon_peri_idx);
X_peri = xxTot(moon_peri_idx,:);
idx_before = find(t_nodes_IG < t_peri);
idx_before = idx_before(end);
t_nodes_IG = [t_nodes_IG(1:idx_before); t_peri; t_nodes_IG(idx_before+1:end)];
X_nodes_IG = [X_nodes_IG(1:idx_before,:); X_peri; X_nodes_IG(idx_before+1:end,:)];
moon_peri_idx_nodes = idx_before + 1;
n_nodes = n_nodes + 1;
if graph
% Show the nodes on the initial guess trajectory
figure;
hold on
box on
grid on
plot3(xxTot(:,1), xxTot(:,2), xxTot(:,3), '-', 'color', 'r','DisplayName', 'Trajectory');
plot3(X_nodes_IG(:,1), X_nodes_IG(:,2), X_nodes_IG(:,3), 'o', 'color','k', 'DisplayName', 'Nodes');
plot3(1-mu_EM, 0, 0, '.', 'color', 'k', 'MarkerSize', 10, 'DisplayName', 'Moon');
plot3(-mu_EM, 0, 0, '.', 'color', 'b', 'MarkerSize', 15, 'DisplayName', 'Earth');
axis equal;
title('Initial Guess Trajectory with Nodes - Adimensional EMB Centered Rotating Frame', 'Interpreter','latex')
set( gca, 'FontName','Times','FontSize',16);
xlabel('x (DU)', 'Interpreter', 'Latex');
ylabel('y (DU)', 'Interpreter', 'Latex');
zlabel('z (DU)', 'Interpreter', 'Latex');
legend('show', 'Interpreter', 'latex');
xlim([min(xxTot(:,1))*1.05 max(xxTot(:,1))*1.5]);
ylim([min(xxTot(:,2)) max(xxTot(:,2))]);
set(gcf,"Theme","light");
end
% Initial guess for regularized impulse variables U = (u, w, s)

% Oshima's regularization:
% Delta v = [dv_x; dv_y; dv_z]
% dv_x = u^2 - w^2 - s^2
% dv_y = 2*u*w
% dv_z = 2*u*s
%
% |Delta v| = u^2 + w^2 + s^2   (objective)
%
% We place impulses at:
%   - Node 1             -> departure from RPO
%   - Node n_nodes       -> perigee of Earth flyby

U_nodes_IG = eps * ones(n_nodes, 3);  % default: almost zero DV everywhere
% 1) Departure DV from RPO 
% Convert DV0 into regularized variables U0 = (u0, w0, s0)
DV0_CR3BP = FullSolution.dv0CR3BP(:);
[U0, is_case2b_0] = dv_to_regularized_U(DV0_CR3BP);

U_nodes_IG(1, :) = U0(:).';  % store as row

% Keep track of nodes that fall in the "case (ii-b)" (if any)
I_case2b = [];
if is_case2b_0
    I_case2b = [I_case2b; 1];
end
% 2) Earth perigee DV 
% Convert DVp into regularized variables Up = (up, wp, sp)
DVp_in = FullSolution.deltavp(:)';
xxp_in = rot2in(X_nodes_IG(end,:), t_nodes_IG(end),mu_EM);
xxp_in = [xxp_in(1:3)*DU, xxp_in(4:6)*VU];
xxpaft_in = [xxp_in(1:3), xxp_in(4:6) + DVp_in];
xxpaft_rot = in2rot([xxpaft_in(1:3)/DU, xxpaft_in(4:6)/VU], t_nodes_IG(end), mu_EM);
DVp_rot = xxpaft_rot(4:6) - X_nodes_IG(end,4:6);
[Up, is_case2b_p] = dv_to_regularized_U(DVp_rot);

U_nodes_IG(end, :) = Up(:).';  % store as row

% Keep track of nodes that fall in the "case (ii-b)" of Oshima (if any)
if is_case2b_p
    I_case2b = [I_case2b; n_nodes];
end

% Pack decision vector y0  (states + regularized impulses + times)
% These are the variable that fmincon will tune

% Layout:
% y = [ x1 y1 z1 vx1 vy1 vz1 u1 w1 s1 ...
%       ...
%       xN yN zN vxN vyN vzN uN wN sN ...
%       t1 ... tN ]   -> total length = 10 * N

y0 = zeros(10 * n_nodes + 1, 1);

% Extract vinf_dir and vInfMag
vInfMag = norm(vInfPlusEarth);
vinf_dir = vInfPlusEarth/vInfMag;

for j = 1:n_nodes
    offset = 9*(j-1);
    y0(offset + 1 : offset + 6) = X_nodes_IG(j, :).';
    y0(offset + 7 : offset + 9) = U_nodes_IG(j, :).';
end

% Append time nodes
y0(9*n_nodes + 1 : 10*n_nodes) = t_nodes_IG;
y0(end) = vInfMag;

% Build constants vector for optimizer
constants = [mu_EM; rminEarth_DU; rminMoon_DU; muE; vinf_dir; xx0RPO; DU; VU; n_nodes; maxTOF; moon_peri_idx_nodes; dv_budget; I_case2b];

% Gradient check
if gradient_check
    opts = optimoptions("fmincon",FiniteDifferenceType="central");
    [validObj, errObj] = checkGradients(@(y) objfun_CR3BP_reachability(y), y0,opts,Display="on");
    
    opts = optimoptions("fmincon",FiniteDifferenceType="central");
    [validCon, errCon] = checkGradients(@(y) constr_CR3BP_reachability(y, constants), y0,opts,IsConstraint = true,Display="on");
    return
end
% fmincon options and call
% Optimization settings
opt = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'MaxIter',  2e3, ...
    'MaxFunEvals', 1e6, ...
    'TolCon', 1e-6, ...
    'TolFun', 1e-6, ...
    'OptimalityTolerance', 1e-2, ...
    'StepTolerance', 1e-4, ...
    'GradObj', 'on', ...      % objective gradient is provided
    'GradConstr', 'on', ...   % constraint gradient is provided
    'DerivativeCheck', 'off', ...
    'FinDiffType','central');

% Call to fmincon:
[y_opt, ~, exitflag1, output1]  = fmincon( ...
    @(y) obj_zero(y), ...
    y0, ...
    [], [], [], [], [], [], ...
    @(y) constr_CR3BP_reachability(y, constants), ...
    opt);

[y_opt, J_opt, exitflag, output2] = fmincon( ...
    @(y) objfun_CR3BP_reachability(y), ...
    y_opt, ...
    [], [], [], [], [], [], ...
    @(y) constr_CR3BP_reachability(y, constants), ...
    opt);

output1.exitflag = exitflag1;
output2.exitflag = exitflag;
output = {output1, output2};

% Plot of the solution
% Store the optimized time nodes
t_nodes = y_opt(9*n_nodes + 1 : 10*n_nodes);
% Extract the optimized state and impulse variables
X_nodes = reshape(y_opt(1:9*n_nodes), [9, n_nodes]);
xx_nodes = X_nodes(1:6,:)';
U_nodes = X_nodes(7:9,:)';
% Extract magnitude of vInf
vInfMag = y_opt(end);

% Compute trajectory in rotating frame
xx_traj_rot = [];
t_traj = [];
STM_vec = cell(n_nodes-1,1);
dimensions = '3D';
framework = 'CR3BP';
stm_computation = 'true'; 
for j = 1:n_nodes-1
    [~,~,STM_hist, x_temp, t_temp] = propagate(xx_nodes(j,:), t_nodes(j), t_nodes(j+1), mu_EM, ...
        dimensions, framework, stm_computation);
    x_temp = x_temp(:,1:6);
    xx_traj_rot = [xx_traj_rot; x_temp];
    t_traj = [t_traj; t_temp];
    STM_vec{j} = STM_hist; 
end

% Transform into ECI
xx_traj_in = rot2in(xx_traj_rot, t_traj, mu_EM);
xx_traj_in(:,1:3) = xx_traj_in(:,1:3) * DU;
xx_traj_in(:,4:6) = xx_traj_in(:,4:6) * VU;

xx_nodes_in = rot2in(xx_nodes, t_nodes, mu_EM);
xx_nodes_in(:,1:3) = xx_nodes_in(:,1:3) * DU;
xx_nodes_in(:,4:6) = xx_nodes_in(:,4:6) * VU;

% Compute impulses
% DV mag
DV_mag_opt = U_nodes(:,1).^2 + U_nodes(:,2).^2 + U_nodes(:,3).^2  ; 
%DV_mag_opt_ms = DV_mag_opt*1000;

% DV vec
DVx_opt = U_nodes(:,1).^2 - U_nodes(:,2).^2 - U_nodes(:,3).^2  ;  
DVy_opt = 2 * U_nodes(:,1) .* U_nodes(:,2)  ;  
DVz_opt = 2 * U_nodes(:,1) .* U_nodes(:,3)  ;  
if ~isempty(I_case2b) 
    DVx_opt(I_case2b) = 2 * U_nodes(I_case2b,1) .* U_nodes(I_case2b,2); 
    DVy_opt(I_case2b) = U_nodes(I_case2b,1).^2 - U_nodes(I_case2b,2).^2 - U_nodes(I_case2b,3).^2;  
    DVz_opt(I_case2b) = 2 * U_nodes(I_case2b,1) .* U_nodes(I_case2b,3);  
end
DVvec_opt_rot = [DVx_opt, DVy_opt, DVz_opt] ;
% DV vec inertial
inertial_no_dv = rot2in(xx_nodes, t_nodes, mu_EM);
inertial_with_dv = rot2in([xx_nodes(:,1:3), xx_nodes(:,4:6) + DVvec_opt_rot], t_nodes, mu_EM);
diff = inertial_with_dv - inertial_no_dv;
DVvec_opt_in = diff(:,4:6) * VU * 1000; % [m/s]

% total DV
dv_opt = (1000*VU)*sum(DV_mag_opt) ;   % [m/s]

% total TOF
tof_opt = (TU/86400)*(t_nodes(end)-t_nodes(1)) ;   % [days]  

% Verify distance from Moon and Earth
[~, ~, is_below_minM] = distance_from_point(xx_traj_rot, 1-mu_EM, rminMoon_DU);
[~, ~, is_below_minE] = distance_from_point(xx_traj_rot, -mu_EM, rminEarth_DU);

% Display results
fprintf('Total Delta V: %.4f m/s\n', dv_opt);
fprintf('Total Time of Flight: %.4f days\n', tof_opt);

if is_below_minM == 0
    fprintf('Moon distance limit respected\n');
else
    fprintf('WARNING: Moon distance limit NOT RESPECTED\n');
end

if is_below_minE == 0
    fprintf('Earth distance limit respected\n');
else
    fprintf('WARNING: Earth distance limit NOT RESPECTED\n');
end
if graph
% Plot the trajectory (rotating frame)
figure;
hold on;
box on;
grid on;
plot3(xx_traj_rot(:,1), xx_traj_rot(:,2), xx_traj_rot(:,3), '-','color','k', 'LineWidth', 1, 'DisplayName', 'Optimized trajectory');
plot3(xx_nodes(:,1), xx_nodes(:,2), xx_nodes(:,3), 'o','color','k', 'DisplayName', 'Nodes');
quiver3(xx_nodes(:,1), xx_nodes(:,2), xx_nodes(:,3), DVvec_opt_rot(:,1), DVvec_opt_rot(:,2), DVvec_opt_rot(:,3), 2, 'Color', 'b', 'AutoScale', 'off', 'DisplayName', '$\Delta$V impulses');
title('Optimized Trajectory - Adimensional EMB Centered Rotating Frame', 'Interpreter', 'Latex')
set(gca, 'FontName','Times','FontSize',15);
xlabel('x (DU)', 'Interpreter', 'Latex');
ylabel('y (DU)', 'Interpreter', 'Latex');
zlabel('z (DU)', 'Interpreter', 'Latex');
axis equal;
xlim([min(xx_traj_rot(:,1)) * 1.1, max(xx_traj_rot(:,1)) * 1.8]);
ylim([min(xx_traj_rot(:,2)) * 1.1, max(xx_traj_rot(:,2)) * 1.1]);
legend('show','Interpreter','latex');
set(gca,'FontSize',15,'FontName','times new roman');
set(gcf,"Theme","light");
end
% Calculate Moon path
thetaM = linspace(0,2*pi,100)';
moonPath = [cos(thetaM), sin(thetaM), zeros(length(thetaM),1)] * DU;

% Plot the trajectory (inertial frame)
figure;
hold on;
box on;
grid on;
plot3(xx_traj_in(:,1), xx_traj_in(:,2), xx_traj_in(:,3), '-','color','k', 'LineWidth', 1, 'DisplayName', 'Optimized trajectory');
plot3(xx_nodes_in(:,1), xx_nodes_in(:,2), xx_nodes_in(:,3), 'o','color','k', 'DisplayName', 'Nodes');
plot3(moonPath(:,1), moonPath(:,2), moonPath(:,3), '--', 'color', [0.5,0.5,0.5], 'DisplayName', 'Moon Path');
quiver3(xx_nodes_in(:,1), xx_nodes_in(:,2), xx_nodes_in(:,3), DVvec_opt_in(:,1), DVvec_opt_in(:,2), DVvec_opt_in(:,3), 2, 'Color', 'b', 'AutoScale', 'off', 'DisplayName', '$\Delta$V impulses');
axis equal;
xlim([min([min(xx_traj_in(:,1)), min(moonPath(:,1))]) * 1.1, max([max(xx_traj_in(:,1)), max(moonPath(:,1))]) * 2.2]);
ylim([min([min(xx_traj_in(:,2)), min(moonPath(:,2))]) * 1.1, max([max(xx_traj_in(:,2)), max(moonPath(:,2))]) * 1.1]);
title('Optimized Trajectory - Earth Centered Inertial Frame', 'Interpreter', 'Latex')
set(gca, 'FontName','Times','FontSize',15);
xlabel('x (km)', 'Interpreter', 'Latex');
ylabel('y (km)', 'Interpreter', 'Latex');
zlabel('z (km)', 'Interpreter', 'Latex');
legend('show','Interpreter','latex');
set(gca,'FontSize',15,'FontName','times new roman');
set(gcf,"Theme","light");
end

% Functions

function [Rp, Vp] = get_perilune_state(v_in, v_out)
%GET_PERILUNE_STATE Recover the lunar periapsis position and velocity from two asymptotes.

% Parameters
muM = 4902.83;
% Turning angle (no need to select for feasible ones because we do that in
% the selection of the initial guess)
dotProd = dot(v_in, v_out);
norm2 = norm(v_in)*norm(v_out);
delta = acos(dotProd/norm2);
% Pericenter radius magnitude
v_mid = (norm(v_in) + norm(v_out)) / 2;
RpMag = muM/v_mid^2 * (1/sin(delta/2) - 1);
% Perilune velocity magnidude (vis-viva equation)
VpMag = sqrt(v_mid^2 + 2*muM/RpMag);
% Position and velocity vectors
Rp = RpMag/(2*v_mid*sin(delta/2)) * (v_in - v_out);
Vp = VpMag/(2*v_mid*cos(delta/2)) * (v_in + v_out);
end

function [value,isterminal,direction] = soi_event_CR3BP(~, x, rSOI, mu)
%SOI_EVENT_CR3BP Stop when the propagated state reaches the prescribed lunar SOI radius.

%SOI_EVENT_CR3BP Stop when |r| hits the SOI radius (either direction).
    posMoon = [1-mu 0 0]';
    value      = norm(x(1:3) - posMoon) - rSOI;
    isterminal = 1;
    direction  = 0; % trigger on both inward/outward crossings
end

function [U, is_case2b] = dv_to_regularized_U(DV)
%DV_TO_REGULARIZED_U  Convert Delta v into regularized variables U = (u,w,s)
% according to Oshima's 3D formulation.
%
% INPUT:
%   DV = [dv_x; dv_y; dv_z]
%
% OUTPUT:
%   U         = [u; w; s]
%   is_case2b = true if this DV falls into Oshima's "case (ii-b)" (indefinite)

dvx = DV(1);
dvy = DV(2);
dvz = DV(3);

dv_mag = sqrt(dvx^2 + dvy^2 + dvz^2);
dv_yz  = dvy^2 + dvz^2;

% Very small threshold to detect zero DV
tol = 1e-14;

is_case2b = false;

if dv_mag < tol
    % Case (iii): zero DV
    fprintf('CASE~3\n');
    u = eps;
    w = eps;
    s = eps;

elseif dvx > 0
    % Case (i): dvx > 0
    fprintf('CASE~1\n');
    u = sqrt(0.5 * (dv_mag + dvx));
    w = dvy / (2*u);
    s = dvz / (2*u);

elseif (dvx <= 0) && (dv_yz > tol)
    % Case (ii-a): dvx <= 0 and (dvy^2 + dvz^2) > 0
    fprintf('CASE~2a\n');
    u = sqrt(0.5 * (dv_mag + dvx));
    w = dvy / (2*u);
    s = dvz / (2*u);

elseif (abs(dvx) > tol) && (dv_yz < tol) && (dvx < 0)
    % Case (ii-b): dvx < 0 and dvy = dvz = 0
    % This is the "indefinite" case that requires special handling in constraints.
    fprintf('CASE~2b\n')
    fprintf('CAUTION!!!\n')
    u = sqrt(-0.5 * dvx);
    w = -sqrt(-0.5 * dvx);
    s = 0.0;
    is_case2b = true;

else
    % Fallback: treat as zero DV (should not happen if cases are exhaustive)
    fprintf('Fallback: CAUTION\n');
    u = eps;
    w = eps;
    s = eps;
end

U = [u; w; s];

end
