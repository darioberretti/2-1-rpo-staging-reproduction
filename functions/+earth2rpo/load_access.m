function access = load_access(repository_root, source)
%LOAD_ACCESS Load reusable Earth-access candidate tables from the selected source.
%
% PROTOTYPE:
%   access = earth2rpo.load_access(repository_root, source)
%
% DESCRIPTION:
%   Resolve earth_access.mat in the reference/generated result branch and
%   load its top-level result with rpo.load_result. Generated input must be
%   computed and saved by the access section before transfer sections use it.
%
% INPUT:
%   repository_root  Repository root path.
%   source           Runner setting zenodo (reference) or generated.
%
% OUTPUT:
%   access           Finalized access result; data.maps.h10000 supplies the
%                      candidate table used for transfer-seed reconstruction.

branch = 'reference';
if source == "generated", branch = 'generated'; end
filename = fullfile(repository_root,'results',branch,'earth_to_rpo','earth_access.mat');
if ~isfile(filename)
    error('earth2rpo:AccessMissing', ...
        ['Earth-access input is missing:\n%s\n', ...
        'For reference input, extract accepted_results.zip into the repository root (see README.md).\n', ...
        'For generated input, compute and save the Earth-access section, ', ...
        'then set experiment_database_source = "generated".'],filename);
end
access = rpo.load_result(filename);
end
