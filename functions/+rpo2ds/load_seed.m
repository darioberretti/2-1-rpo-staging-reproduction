function seed = load_seed(repository_root)
%LOAD_SEED Load the tracked planar RPO seed used by Run 02.
%
% PROTOTYPE:
%   seed = rpo2ds.load_seed(repository_root)
%
% DESCRIPTION:
%   Read data/inputs/reference_rpo_seed.mat and require state_2d in the
%   reference_rpo_seed variable. Preserve the seed and its source metadata.
%
% INPUT:
%   repository_root  Repository root path.
%
% OUTPUT:
%   seed             Stored struct with state_2d = [x; y; vx; vy], using
%                      normalized rotating-frame DU and VU.

filename = fullfile(repository_root, 'data', 'inputs', 'reference_rpo_seed.mat');
if ~isfile(filename)
    error('rpo2ds:SeedMissing', [ ...
        'The tracked RPO seed is missing:\n%s\n\n', ...
        'Restore data/inputs/reference_rpo_seed.mat from the repository.'], filename);
end
loaded = load(filename, 'reference_rpo_seed');
if ~isfield(loaded, 'reference_rpo_seed') || ...
        ~isfield(loaded.reference_rpo_seed, 'state_2d')
    error('rpo2ds:SeedSchema', ...
        'The RPO seed MAT has an incompatible schema:\n%s', filename);
end
seed = loaded.reference_rpo_seed;
end
