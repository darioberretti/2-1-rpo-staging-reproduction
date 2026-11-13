function result = compute_access(cfg, inputs)
%COMPUTE_ACCESS Search the historical outer grid and retain inner-boundary hits.
%
% PROTOTYPE:
%   result = earth2rpo.compute_access(cfg, inputs)
%
% DESCRIPTION:
%   Sample the impulse magnitude/direction grid, propagate backward to the
%   ordered boundary altitudes, and pass only successful candidates to each
%   inner search. Preserve the original impulse conversion operation order
%   and unscaled historical C3-like quantity. Repropagate the frozen examples.
%
% INPUT:
%   cfg     cfg.earth_to_rpo; uses shared normalization and access grid,
%             horizon, boundary altitudes, clearances, and parallel settings.
%   inputs  Tracked Run 03 inputs with x0_rpo_bcr4bp [4 x 1] in DU/VU and
%             example_grid_indices identifying the manuscript connections.
%
% OUTPUT:
%   result  Struct with data.maps.h100000/h10000/h1000 candidate tables,
%             data.examples, and data.rpo_reference. Table time_tu is signed
%             backward time [TU]; delta_v_tu retains the historical field name
%             for normalized impulse [VU]. c3_historical is not scaled to
%             km^2/s^2. metrics.candidate_counts follows outer-to-inner order.

% Choose serial/parallel execution before allocating the historical impulse grid.
parallel = rpo.parallel_available(cfg.access.use_parallel);
if cfg.access.use_parallel && ~parallel
    fprintf('Parallel Computing Toolbox is unavailable: using slower serial execution. Install the toolbox to enable parallel computation.\n');
end
n = cfg.access.grid_size;
% Preserve the historical conversion's operation order, not just its units.
maximum = cfg.access.maximum_impulse_mps/1000/cfg.shared.DU*3600*24*cfg.shared.TU_days;
magnitudes = linspace(0,maximum,n);
angles = linspace(0,2*pi,n);
[magnitude_index,angle_index] = ndgrid(1:n,1:n);
map = table(magnitude_index(:),angle_index(:),magnitudes(magnitude_index(:)).', ...
    angles(angle_index(:)).','VariableNames', ...
    {'magnitude_index','angle_index','delta_v_tu','angle_rad'});

% Pass only successful outer-boundary candidates to each subsequent inner search.
for altitude = cfg.access.boundary_altitudes_km
    count = height(map);
    values = NaN(count,2);
    started = tic;
    % Small blocks bound memory and provide an early throughput observation.
    for first = 1:1000:count
        last = min(first+999,count);
        dv = map.delta_v_tu(first:last); theta = map.angle_rad(first:last);
        block = NaN(numel(dv),2);
        if parallel
            parfor k = 1:numel(dv)
                block(k,:) = evaluate(dv(k),theta(k),altitude,cfg,inputs.x0_rpo_bcr4bp);
            end
        else
            for k = 1:numel(dv)
                block(k,:) = evaluate(dv(k),theta(k),altitude,cfg,inputs.x0_rpo_bcr4bp);
            end
        end
        values(first:last,:) = block;
        if first == 1 || last == count || floor(last/(.05*count)) > floor((first-1)/(.05*count))
            fprintf('Earth boundary %.0f km: %.1f%% (%d/%d), elapsed %.1f s\n', ...
                altitude,100*last/count,last,count,toc(started));
        end
    end

    % Store signed hit times and retain finite boundary crossings without reordering.
    map.time_tu = values(:,1); map.c3_historical = values(:,2);
    map = map(isfinite(map.time_tu),:);
    result.data.maps.(sprintf('h%d',altitude)) = map;
end

% Recover the fixed manuscript examples and a reference orbit for rendering.
result.data.examples = earth2rpo.make_examples(map,cfg,inputs);
[~,~,~,result.data.rpo_reference] = propagate(inputs.x0_rpo_bcr4bp,0,10, ...
    cfg.shared.mu,'2D','BCR4BP','false',false,[],pi);
result.metrics.candidate_counts = [height(result.data.maps.h100000), ...
    height(result.data.maps.h10000),height(result.data.maps.h1000)];
end

function values = evaluate(magnitude,angle,altitude,cfg,rpo_state)
%EVALUATE Measure one accepted hit's signed time and historical C3.
%
% DESCRIPTION:
%   Returns [time, c3], or NaNs on a miss. time is backward TU; c3 uses
% normalized Earth-inertial velocity and is not scaled to km^2/s^2.

[hit,time,state] = earth2rpo.propagate_boundary(magnitude,angle,altitude,cfg,rpo_state);
values = [NaN,NaN];
if hit
    mu = cfg.shared.mu;
    radius = sqrt((state(1)+mu)^2+state(2)^2);
    c3 = (state(3)-state(2))^2+(state(4)+state(1)+mu)^2-2*(1-mu)/radius;
    values = [time,c3];
end
end
