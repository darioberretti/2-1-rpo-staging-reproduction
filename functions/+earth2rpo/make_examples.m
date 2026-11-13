function examples = make_examples(map, cfg, inputs)
%MAKE_EXAMPLES Repropagate the frozen manuscript Earth-access connections.
%
% PROTOTYPE:
%   examples = earth2rpo.make_examples(map, cfg, inputs)
%
% DESCRIPTION:
%   Locate the stored magnitude/angle grid pairs in the retained 1000 km
%   boundary map and recover their full backward trajectories.
%
% INPUT:
%   map     Retained inner-boundary candidate table from compute_access.
%   cfg     cfg.earth_to_rpo; uses shared scales and propagation settings.
%   inputs  Tracked x0_rpo_bcr4bp [DU/VU] and example_grid_indices.
%
% OUTPUT:
%   examples  Struct array with grid indices, backward time vectors [TU],
%               rotating state rows [DU/VU], positive TOF [days], historical
%               unscaled C3, and impulse magnitude [m/s].

for k = 1:size(inputs.example_grid_indices,1)
    key = inputs.example_grid_indices(k,:);
    row = find(map.magnitude_index == key(1) & map.angle_index == key(2));
    [~,~,~,time,state] = earth2rpo.propagate_boundary(map.delta_v_tu(row), ...
        map.angle_rad(row), 1000, cfg, inputs.x0_rpo_bcr4bp);
    examples(k) = struct('grid_indices', key, 'time_tu', time, 'state_rot', state, ...
        'tof_days', -map.time_tu(row)*cfg.shared.TU_days, ...
        'c3_historical', map.c3_historical(row), ...
        'impulse_mps', map.delta_v_tu(row)/cfg.shared.TU_days/3600/24*cfg.shared.DU*1000); %#ok<AGROW>
end
end
