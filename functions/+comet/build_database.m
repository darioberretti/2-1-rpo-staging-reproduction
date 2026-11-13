function database = build_database(settings, constants)
%BUILD_DATABASE Compute Earth, L1 and L2 cumulative heliocentric cost grids.
%
% PROTOTYPE:
%   database = comet.build_database(settings, constants)
%
% DESCRIPTION:
%   Sweep intercept radius, rotating angle and available flight time.
%   Departures are fixed rotating locations on the Sun-Earth line.
%   Angle computations use thread workers when requested and available.
%
% INPUT:
%   settings   cfg.comet_case_study.heliocentric: Rc_grid_au [AU],
%              theta_grid_deg [deg], tof_grid_days [days], max_revolutions
%              and use_parallel.
%   constants  cfg.comet_case_study.constants: AU_km, mu_sun_km3s2,
%              day_s, year_days, year_s, n_earth_rads [rad/s], and
%              se_libration_offset_au [AU].
%
% OUTPUT:
%   database   Struct containing the three shared axis vectors and
%              cumulative_best_vinf_mag_kms [R x A x T x 3 double, km/s].
%              The final axis follows departure_modes = Earth, L1, L2.
%              departure_rot_au [3 x 3] stores rotating positions [AU]
%              in columns. Metadata records the model, units, constants,
%              max_revolutions and schema_version. No MAT file is saved here.

timer = tic;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'third_party', 'FEX-Lambert'));
radius = settings.Rc_grid_au(:);
angles = settings.theta_grid_deg(:);
times = settings.tof_grid_days(:);
shape = [numel(radius), numel(angles), numel(times)];
offset = constants.se_libration_offset_au;
modes = ["Earth", "L1", "L2"];

% Keep interpolation axes and generation metadata with the cost tensor.
database = struct('schema_version', 1, 'geometry_model', "rotating_frame_boundaries", ...
    'departure_modes', modes, ...
    'units', struct('radius', 'AU', 'angle', 'degrees', 'tof', 'days', 'excess_speed', 'km/s'));
database.Rc_grid_au = radius;
database.theta_grid_deg = angles;
database.tof_grid_days = times;
database.cumulative_best_vinf_mag_kms = NaN([shape, 3]);
database.departure_rot_au = [1, 1-offset, 1+offset; zeros(2, 3)];
database.max_revolutions = settings.max_revolutions;
for name = {'AU_km','mu_sun_km3s2','day_s','year_days','year_s', ...
        'n_earth_rads','se_libration_offset_au'}
    database.constants.(name{1}) = constants.(name{1});
end

% Prepare workers before evaluating independent angle pairs.
parallel = rpo.parallel_available(settings.use_parallel);
if settings.use_parallel && ~parallel
    fprintf('Parallel Computing Toolbox is unavailable; rebuilding serially will take longer. Install it for faster rebuilding.\n');
end
if parallel && isempty(gcp('nocreate')), parpool('Threads'); end

for origin = 1:3
    settings.departure_rot_au = database.departure_rot_au(:, origin);
    for i_radius = 1:shape(1)
        pairs = cell(shape(2), 1);
        radius_value = radius(i_radius);
        if parallel
            parfor i_angle = 1:shape(2)
                pairs{i_angle} = comet.build_pair(radius_value, angles(i_angle), settings, constants);
            end
        else
            for i_angle = 1:shape(2)
                pairs{i_angle} = comet.build_pair(radius_value, angles(i_angle), settings, constants);
            end
        end

        for i_angle = 1:shape(2)
            database.cumulative_best_vinf_mag_kms(i_radius, i_angle, :, origin) = pairs{i_angle};
        end
        fprintf('%s Lambert database: %d/%d radius slices (%.1f%%), %.1f min elapsed.\n', ...
            modes(origin), i_radius, shape(1), 100*i_radius/shape(1), toc(timer)/60);
    end
end
end
