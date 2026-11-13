function [inputs, metadata] = load_reachability_initial_guesses(repository_root)
%LOAD_REACHABILITY_INITIAL_GUESSES Load the frozen historical preprocessing.
%
% PROTOTYPE:
%   [inputs, metadata] = rpo2ds.load_reachability_initial_guesses(repository_root)
%
% DESCRIPTION:
%   Read the tracked run_02_reachability_initial_guesses.mat input. This
%   supplies the historical patched starts; the directional optimizations
%   still execute. compute_reachability validates their grid/configuration.
%
% INPUT:
%   repository_root  Repository root path.
%
% OUTPUT:
%   inputs           Stored reachability_initial_guesses struct with
%                      schema_version, configuration, provenance, rpo_state_2d,
%                      and ordered points containing the patched starts.
%   metadata         Logical name, frozen-input label, schema version, and hash.

filename = fullfile(repository_root, 'data', 'inputs', ...
    'run_02_reachability_initial_guesses.mat');
if ~isfile(filename)
    error('rpo2ds:ReachabilityInitialGuessesMissing', ...
        ['The tracked Run 02 reachability preprocessing is missing:\n%s\n\n', ...
        'Restore data/inputs/run_02_reachability_initial_guesses.mat ', ...
        'from the repository.'], filename);
end
loaded = load(filename, 'reachability_initial_guesses');
if ~isfield(loaded, 'reachability_initial_guesses')
    error('rpo2ds:ReachabilityInitialGuessesSchema', ...
        'The Run 02 reachability preprocessing has an incompatible schema:\n%s', ...
        filename);
end
inputs = loaded.reachability_initial_guesses;
metadata = struct('logical_name', "run_02_reachability_initial_guesses.mat", ...
    'artifact_type', "frozen_input", 'schema_version', inputs.schema_version, ...
    'sha256', rpo2ds.sha256_file(filename));
end
