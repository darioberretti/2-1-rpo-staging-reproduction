function result = finalize_result(result, cfg, producer, artifact_type, input_metadata)
%FINALIZE_RESULT Add the compact Run 02 provenance and completion status.
%
% PROTOTYPE:
%   result = rpo2ds.finalize_result(result, cfg, producer, artifact_type, input_metadata)
%
% DESCRIPTION:
%   Attach the Run 02 configuration, UTC creation time, Git revision, MATLAB
%   and toolbox versions, and supplied input metadata. Preserve existing
%   study-specific provenance and status. Completion is not an acceptance test.
%
% INPUT:
%   result          Compact study result struct.
%   cfg             Full reproduction_config() output.
%   producer        Producing runner/section label.
%   artifact_type   reference or generated.
%   input_metadata  Struct of source inputs, database roles, and hashes.
%
% OUTPUT:
%   result          Annotated struct with configuration = cfg.rpo_to_deep_space,
%                   shared provenance, and completion/failure-status fields.

% Capture the scientific inputs, environment, and current source revision.
result.configuration = cfg.rpo_to_deep_space;
if ~isfield(result, 'provenance') || ~isstruct(result.provenance)
    result.provenance = struct();
end
result.provenance.schema_version = "1.0";
result.provenance.artifact_type = string(artifact_type);
result.provenance.producer = string(producer);
result.provenance.created_utc = string(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
result.provenance.code_version = local_code_version();
result.provenance.matlab_version = string(version);
result.provenance.toolbox_versions = local_toolbox_versions();
result.provenance.input_metadata = input_metadata;

% Preserve producer status, supplying only absent completion/failure fields.
if ~isfield(result, 'status') || ~isstruct(result.status)
    result.status = struct('completed', true, 'failed_sections', strings(0, 1));
elseif ~isfield(result.status, 'completed')
    result.status.completed = true;
end
if ~isfield(result.status, 'failed_sections')
    result.status.failed_sections = strings(0, 1);
end
end

function versions = local_toolbox_versions()
%LOCAL_TOOLBOX_VERSIONS Record installed optimization/parallel versions and CasADi presence.

versions = struct();
for product = ["optim", "parallel"]
    information = ver(char(product));
    field = char(product);
    if isempty(information)
        versions.(field) = "not installed";
    else
        versions.(field) = string(information(1).Version);
    end
end
versions.casadi_ipopt_available = ~isempty(which('casadi.Opti'));
end

function code_version = local_code_version()
%LOCAL_CODE_VERSION Read HEAD and append +uncommitted when the working-tree query reports changes.

[status, revision] = system('git rev-parse --short HEAD');
if status ~= 0
    code_version = "unavailable";
    return
end
code_version = string(strtrim(revision));
[dirty_status, changes] = system('git status --porcelain --untracked-files=all');
if dirty_status == 0 && strlength(strtrim(string(changes))) > 0
    code_version = code_version + "+uncommitted";
end
end
