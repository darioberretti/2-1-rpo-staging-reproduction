function targets = load_targets(repository_root)
%LOAD_TARGETS Load the compact accepted Venus/Mars target inputs.
%
% PROTOTYPE:
%   targets = rpo2ds.load_targets(repository_root)
%
% DESCRIPTION:
%   Read data/inputs/rpo_to_deep_space_targets.mat and require venus/mars
%   structs. The consumers use stored planetary states and Earth-relative
%   target velocities; no ephemeris query or target redesign occurs here.
%
% INPUT:
%   repository_root  Repository root path.
%
% OUTPUT:
%   targets          Stored struct. Each planet supplies departure/arrival
%                      epochs [MJD2000 days], flight time [days], heliocentric
%                      positions [km], Lambert departure velocity [km/s], and
%                      Earth-relative v-infinity in ECI/RTH [km/s]. The parent
%                      struct supplies the solar gravitational parameter.

filename = fullfile(repository_root, 'data', 'inputs', ...
    'rpo_to_deep_space_targets.mat');
if ~isfile(filename)
    error('rpo2ds:TargetsMissing', [ ...
        'The compact target input is missing:\n%s\n\n', ...
        'Restore data/inputs/rpo_to_deep_space_targets.mat from the repository.'], ...
        filename);
end
loaded = load(filename, 'targets');
if ~isfield(loaded, 'targets') || ~isstruct(loaded.targets) || ...
        ~all(isfield(loaded.targets, {'venus', 'mars'}))
    error('rpo2ds:TargetsSchema', ...
        'The compact target input has an incompatible schema:\n%s', filename);
end
targets = loaded.targets;
end
