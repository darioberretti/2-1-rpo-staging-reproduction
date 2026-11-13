function result = compute_accessibility(settings,directory)
%COMPUTE_ACCESSIBILITY Extract full-range and three limited-time heliocentric cost maps.
%
% PROTOTYPE:
%   result = comet.compute_accessibility(settings, directory)
%
% DESCRIPTION:
%   Read cumulative minimum excess-speed grids for Earth, L1,
%   L2 and their envelope. The full-range slice uses the final stored
%   time, while limited horizons must match stored flight-time entries.
%
% INPUT:
%   settings  Comet-study cfg.comet_case_study struct returned by
%               reproduction_config(); see config/reproduction_config.m.
%               Required substructures are identified below.
%   settings.heliocentric.limited_tof_days lists three horizons [days].
%   directory  Folder containing the comet.database_file MAT product.
%
% OUTPUT:
%   result    Struct with settings, units and maps [1 x 4] in Earth,
%               L1, L2, L1/L2 order. Each map contains label, radius_au
%               [R x 1, AU], angle_deg [A x 1, deg], horizons_days
%               [1 x 4, days], cost_kms [R x A x 4, km/s] and database
%               identity (filename/bytes/sha256). Slice 1 is full range.

result = struct('settings',settings.heliocentric,'maps',struct([]));
maps = cell(1,4);
labels = ["Earth","L1","L2","L1/L2"];
for k = 1:4
    database = comet.read_database(directory,k);
    horizons = [database.tof_grid_days(end),settings.heliocentric.limited_tof_days];
    indices = arrayfun(@(time) find(database.tof_grid_days==time,1),horizons);
    maps{k} = struct('label',labels(k),'radius_au',database.Rc_grid_au, ...
        'angle_deg',database.theta_grid_deg,'horizons_days',horizons, ...
        'cost_kms',database.cumulative_best_vinf_mag_kms(:,:,indices), ...
        'database',database.identity);
    fprintf('Accessibility slices loaded: %s, full range %.0f days.\n',labels(k),horizons(1));
end
result.maps = [maps{:}];
result.units = struct('radius','AU','angle','degrees','tof','days','excess_speed','km/s');
end
