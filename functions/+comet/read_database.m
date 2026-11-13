function database = read_database(directory, index)
%READ_DATABASE Read one departure role's cumulative heliocentric cost grid.
%
% PROTOTYPE:
%   database = comet.read_database(directory, index)
%
% DESCRIPTION:
%   Read the shared axes and the requested Earth, L1 or L2 cost slice.
%   For the combined role, select L2 when L1 is nonfinite or L2 is finite
%   and strictly cheaper. Record the file identity and departure role.
%
% INPUT:
%   directory  Folder containing the comet.database_file MAT product.
%   index      Role index: 1 Earth, 2 L1, 3 L2 or 4 combined L1/L2.
%
% OUTPUT:
%   database   Struct with Rc_grid_au [R x 1, AU], theta_grid_deg
%                [A x 1, deg], tof_grid_days [T x 1, days], and
%                cumulative_best_vinf_mag_kms [R x A x T, km/s]. identity
%                records filename, bytes, sha256 and departure_mode.

filename = comet.database_file();
path = fullfile(directory, filename);
if ~isfile(path)
    error('comet:DatabaseMissing',['Heliocentric database is missing:\n%s\n', ...
        'Extract comet_lambert_databases.zip for reference data (see docs/DATA.md), ', ...
        'or run the Run 06 accessibility section with heliocentric_database_mode = "rebuild".'],path);
end

% Read the common interpolation axes and one origin slice.
database = struct();
for name = ["Rc_grid_au", "theta_grid_deg", "tof_grid_days"]
    database.(name) = h5read(path, "/database/" + name);
end
dataset = '/database/cumulative_best_vinf_mag_kms';
info = h5info(path, dataset);
count = [info.Dataspace.Size(1:3), 1];
origin = index;
if index == 4, origin = 2; end
costs = h5read(path, dataset, [1, 1, 1, origin], count);

% Form the combined L1/L2 envelope.
if index == 4
    l2 = h5read(path, dataset, [1, 1, 1, 3], count);
    use_l2 = ~isfinite(costs) | (isfinite(l2) & l2 < costs);
    costs(use_l2) = l2(use_l2);
end
database.cumulative_best_vinf_mag_kms = costs;
roles = ["Earth", "L1", "L2", "L1/L2"];
database.identity = struct('filename', filename, 'bytes', dir(path).bytes, ...
    'sha256', rpo2ds.sha256_file(path), 'departure_mode', roles(index));
end
