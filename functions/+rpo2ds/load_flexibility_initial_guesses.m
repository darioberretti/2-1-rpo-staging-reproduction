function [preprocessing, metadata] = load_flexibility_initial_guesses(repository_root)
%LOAD_FLEXIBILITY_INITIAL_GUESSES Load archived Run 02 optimizer setups.
%
% PROTOTYPE:
%   [preprocessing, metadata] = rpo2ds.load_flexibility_initial_guesses(repository_root)
%
% DESCRIPTION:
%   Read the tracked run_02_flexibility_initial_guesses.mat input and check
%   its container fields. compute_flexibility checks the grid/recipe match;
%   solve_frozen_flexibility_case performs the subsequent optimization.
%
% INPUT:
%   repository_root  Repository root path.
%
% OUTPUT:
%   preprocessing    Stored struct with schema_version, configuration, provenance,
%                      and venus/mars untuned_setup and retuning_setup cells.
%                      Each setup contains y0, lb, ub, problem_data, and
%                      velocity/time scales.
%   metadata         Logical name, tracked-preprocessing label, and SHA-256.

filename = fullfile(repository_root, 'data', 'inputs', ...
    'run_02_flexibility_initial_guesses.mat');
if ~isfile(filename)
    error('rpo2ds:FlexibilityPreprocessingMissing', [ ...
        'The tracked Run 02 flexibility preprocessing is missing:\n%s\n\n', ...
        'Restore data/inputs/run_02_flexibility_initial_guesses.mat ', ...
        'from the repository.'], filename);
end
loaded = load(filename, 'preprocessing');
if ~isfield(loaded, 'preprocessing') || ...
        ~all(isfield(loaded.preprocessing, {'schema_version', 'configuration', ...
        'provenance', 'venus', 'mars'}))
    error('rpo2ds:FlexibilityPreprocessingSchema', ...
        'The Run 02 flexibility preprocessing has an incompatible schema:\n%s', ...
        filename);
end
preprocessing = loaded.preprocessing;
metadata = struct('logical_name', "run_02_flexibility_initial_guesses.mat", ...
    'artifact_type', "tracked preprocessing", ...
    'sha256', rpo2ds.sha256_file(filename));
end
