function seed = build_seed(map, cfg, inputs)
%BUILD_SEED Reconstruct and rank the historical Earth-side conic patches.
%
% PROTOTYPE:
%   seed = earth2rpo.build_seed(map, cfg, inputs)
%
% DESCRIPTION:
%   Repropagate candidate boundary states, recover Earth-centered two-body
%   periapses, and retain outward/positive-time branches in the prescribed
%   altitude range. Rank by circular-speed mismatch, historical C3, then
%   flight time. Join the selected conic and reversed database arc as a seed;
%   this is the active historical path without the extra Pareto experiments.
%
% INPUT:
%   map     10000 km candidate table from the selected saved access result.
%   cfg     cfg.earth_to_rpo with shared, access, patch, and validation settings.
%   inputs  Tracked x0_rpo_bcr4bp [4 x 1] in rotating DU/VU.
%
% OUTPUT:
%   seed    Struct with constants, selected_candidate, match/periapsis states,
%             conic/database/concatenated timelines and rotating trajectories,
%             plot_data, and summary_counts. Rotating states use DU/VU and
%             times TU; inertial states use km/km/s. Selected-candidate fields
%             report periapsis altitude [km], mismatch [km/s], TOF [days],
%             and unscaled historical C3.

% Allocate one reconstructed conic state and ranking record per saved candidate.
s = cfg.shared;
p = cfg.patch;
n = height(map);
match_time = NaN(n,1); match_state = NaN(n,4); match_inertial = NaN(n,4);
periapsis_altitude = NaN(n,1); mismatch = NaN(n,1);
positive_branch = false(n,1);

% Recover each matching-boundary state and its outward Earth-periapsis conic.
for k = 1:n
    [hit,time,state] = earth2rpo.propagate_boundary(map.delta_v_tu(k), ...
        map.angle_rad(k),p.matching_altitude_km,cfg,inputs.x0_rpo_bcr4bp);
    if hit
        inertial = rot2in(state,time,s.mu);
        inertial(1:2) = inertial(1:2)*s.DU;
        inertial(3:4) = inertial(3:4)*s.VU;
        [rp,vp,conic] = planar_conic_pericenter_from_match_state( ...
            inertial(1:2),inertial(3:4),398600,0);
        radial_direction = rp/norm(rp);
        tangent = sign(conic.h)*[-radial_direction(2),radial_direction(1)];
        circular_velocity = sqrt(398600/conic.rp_radius)*tangent;
        match_time(k) = time; match_state(k,:) = state;
        match_inertial(k,:) = inertial;
        periapsis_altitude(k) = conic.rp_altitude;
        mismatch(k) = norm(vp-circular_velocity);
        positive_branch(k) = conic.tof_peri_to_match > 0 && conic.radial_velocity_match > 0;
    end
    if mod(k,500) == 0 || k == n
        fprintf('Earth conic reconstruction: %d/%d candidates\n',k,n);
    end
end

% Filter the physical branch and altitude range, then apply the historical ranking.
valid = positive_branch & periapsis_altitude >= p.periapsis_altitude_range_km(1) ...
    & periapsis_altitude <= p.periapsis_altitude_range_km(2);
indices = find(valid);
ranking = [mismatch(valid),map.c3_historical(valid),-map.time_tu(valid)*s.TU_days];
[~,order] = sortrows(ranking,[1,2,3]);
selected = indices(order(1));

% Reconstruct the selected conic arc and align it with the reversed database leg.
[rp,vp,conic] = planar_conic_pericenter_from_match_state( ...
    match_inertial(selected,1:2),match_inertial(selected,3:4),398600,p.conic_samples);
[~,~,~,database_time,database_state] = earth2rpo.propagate_boundary( ...
    map.delta_v_tu(selected),map.angle_rad(selected),p.matching_altitude_km,cfg,inputs.x0_rpo_bcr4bp);
conic_start_time = match_time(selected) - conic.tof_peri_to_match/(s.TU_days*86400);
conic_time = conic_start_time + conic.t_arc/(s.TU_days*86400);
conic_state = in2rot([conic.rr_arc/s.DU,conic.vv_arc/s.VU],conic_time,s.mu);
database_time = flipud(database_time); database_state = flipud(database_state);
[~,~,~,rpo_state] = propagate(inputs.x0_rpo_bcr4bp,0,s.T_RPO,s.mu, ...
    '2D','BCR4BP','false',false,0,0);

% Prepare the departure-circle context and package reusable seed/plot data.
circle_angle = linspace(0,2*pi,400);
circle_radius = (s.R_e+cfg.validation.departure_altitude_km)/s.DU;
circle_x = -s.mu + circle_radius*cos(circle_angle);
circle_y = circle_radius*sin(circle_angle);
seed.constants = s;
seed.selected_candidate = table(selected,map.delta_v_tu(selected),map.angle_rad(selected),0, ...
    map.c3_historical(selected),-map.time_tu(selected)*s.TU_days, ...
    periapsis_altitude(selected),mismatch(selected), ...
    'VariableNames',{'source_index','deltaV','thetaV','thetaSun0','C3','TOF_days', ...
    'periapsis_altitude_km','vp_circular_mismatch_kms'});
seed.match_time_tu = match_time(selected);
seed.match_state_rot = match_state(selected,:);
seed.match_state_inertial = match_inertial(selected,:);
seed.periapsis_state_inertial = [rp,vp];
seed.database_time_tu = database_time; seed.database_state_rot = database_state;
seed.conic_time_tu = conic_time; seed.conic_state_rot = conic_state;
seed.concatenated_time_tu = [conic_time(:);database_time(2:end)];
seed.concatenated_state_rot = [conic_state;database_state(2:end,:)];
seed.plot_data = struct('mu',s.mu,'xx_rpo',rpo_state,'vleo_x',circle_x, ...
    'vleo_y',circle_y,'xx_db',database_state,'xx_conic',conic_state,'moon_xy',[1-s.mu,0]);
seed.summary_counts = struct('total_candidates',n,'positive_branch_candidates',nnz(positive_branch), ...
    'altitude_range_candidates',nnz(valid));
fprintf('Selected patch: source row %d, periapsis %.9f km, departure mismatch %.6f m/s\n', ...
    selected,periapsis_altitude(selected),1000*mismatch(selected));
end
