function [database, metadata] = load_database(repository_root, source, database_role)
%LOAD_DATABASE Load one reference/generated RPO-to-Moon database.
%
% PROTOTYPE:
%   [database, metadata] = rpo2ds.load_database(repository_root, source, database_role)
%
% DESCRIPTION:
%   Resolve the table file by source and scientific role, require one table
%   variable and the core schema columns, then record its source and hash.
%   Loading does not convert either historical database into the other schema.
%   The two table roles are described in docs/DATA.md.
%
% INPUT:
%   repository_root  Repository root path.
%   source           zenodo (data/reference) or generated (data/generated).
%   database_role    phase_aware or legacy_dense.
%
% OUTPUT:
%   database         Table containing vInfMagMoon [km/s], vInfMoon_RTH, deltaV,
%                      and t_event_days [days]. Phase-aware tables also require
%                      orbitPhase [TU] and n_hits. Remaining consumer columns
%                      retain their source schema and units.
%   metadata         Logical name, source/type, schema, row count, columns,
%                      and SHA-256 digest.

% Resolve the requested source and scientific database role.
source = validatestring(char(source), {'zenodo', 'generated'});
database_role = validatestring(char(database_role), {'phase_aware', 'legacy_dense'});
if strcmp(source, "zenodo")
    branch = 'reference';
else
    branch = 'generated';
end
if strcmp(database_role, "phase_aware")
    filename = 'RPOMoonDatabase10x3_50phases.mat';
else
    filename = 'RPOMoonDatabase50x3.mat';
end
path = fullfile(repository_root, 'data', branch, 'rpo_to_moon', filename);
if ~isfile(path)
    if strcmp(source, "zenodo")
        recovery = 'Extract the relevant database ZIP into the repository root (see docs/DATA.md).';
    else
        recovery = 'Run data_builders/build_rpo_to_moon_database.m and select the generated database.';
    end
    error('rpo2ds:DatabaseMissing', ...
        'The %s %s RPO-to-Moon database is missing:\n%s\n\n%s', ...
        source, database_role, path, recovery);
end

% Require one table and the core columns needed by the selected role.
loaded = load(path);
names = fieldnames(loaded);
table_names = names(cellfun(@(name) istable(loaded.(name)), names));
if numel(table_names) ~= 1
    error('rpo2ds:DatabaseSchema', ...
        'The selected RPO-to-Moon database must contain one MATLAB table:\n%s', path);
end
database = loaded.(table_names{1});
required = {'vInfMagMoon', 'vInfMoon_RTH', 'deltaV', 't_event_days'};
if strcmp(database_role, 'phase_aware')
    required = [required, {'orbitPhase', 'n_hits'}];
end
missing = required(~ismember(required, database.Properties.VariableNames));
if ~isempty(missing)
    error('rpo2ds:DatabaseSchema', ...
        'The selected %s RPO-to-Moon database is missing columns %s:\n%s', ...
        database_role, strjoin(missing, ', '), path);
end

% Record the source/schema without replacing historical normalization or phase data.
metadata = struct('logical_name', string(filename), 'source', string(source), ...
    'artifact_type', string(local_artifact_type(source)), ...
    'schema', string(database_role), ...
    'row_count', height(database), ...
    'variable_names', string(database.Properties.VariableNames), ...
    'sha256', rpo2ds.sha256_file(path));
end

function artifact_type = local_artifact_type(source)
%LOCAL_ARTIFACT_TYPE Map the accepted/generated database source to its provenance label.

if strcmp(source, 'zenodo')
    artifact_type = "reference";
else
    artifact_type = "generated";
end
end
