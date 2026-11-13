function metadata = input_metadata(repository_root, include_targets)
%INPUT_METADATA Return logical names and hashes for tracked Run 02 inputs.
%
% PROTOTYPE:
%   metadata = rpo2ds.input_metadata(repository_root, include_targets)
%
% DESCRIPTION:
%   Record the reference RPO seed and, when requested, the accepted target
%   file without copying their payload into the result provenance.
%
% INPUT:
%   repository_root  Repository root path.
%   include_targets  Optional logical, default true.
%
% OUTPUT:
%   metadata         Struct with seed and optional targets entries, each
%                      containing logical_name and the SHA-256 digest or the
%                      missing/unreadable sentinel from rpo2ds.sha256_file.

if nargin < 2
    include_targets = true;
end
metadata.seed = local_one(repository_root, 'data/inputs/reference_rpo_seed.mat');
if include_targets
    metadata.targets = local_one(repository_root, 'data/inputs/rpo_to_deep_space_targets.mat');
end
end

function item = local_one(repository_root, relative_name)
%LOCAL_ONE Pair one repository-relative input name with its file digest.

filename = fullfile(repository_root, strrep(relative_name, '/', filesep));
item = struct('logical_name', string(relative_name), ...
    'sha256', rpo2ds.sha256_file(filename));
end
