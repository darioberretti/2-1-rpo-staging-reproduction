function result = finalize_result(result, cfg, producer, input_sources, artifact_type)
%FINALIZE_RESULT Add shared configuration, provenance, and status fields.
%
% PROTOTYPE:
%   result = rpo.finalize_result(result, cfg, producer, input_sources, artifact_type)
%
% DESCRIPTION:
%   Attach the scientific configuration, UTC creation time, MATLAB release,
%   Git revision, and schema version. Mark completion and artifact type.
%   The completion flag records producer completion, not scientific validation.
%
% INPUT:
%   result         Study result struct to annotate.
%   cfg            Configuration from reproduction_config(), including meta.name.
%   producer       Producing script or function name.
%   input_sources  Text or text array identifying source inputs/settings.
%   artifact_type  Artifact label, normally reference or generated.
%
% OUTPUT:
%   result         Input struct with configuration, provenance, and status.

% Capture the scientific inputs and the environment that produced the result.
result.configuration = cfg;
result.provenance.producer = string(producer);
result.provenance.configuration_name = string(cfg.meta.name);
result.provenance.input_sources = string(input_sources);
result.provenance.created_utc = string(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
result.provenance.matlab_release = string(version('-release'));
result.provenance.code_version = local_code_version();
result.provenance.schema_version = "1.0";

% Record completion separately from any numerical comparison or acceptance test.
result.status.completed = true;
result.status.artifact_type = string(artifact_type);
end

function code_version = local_code_version()
%LOCAL_CODE_VERSION Read the current Git revision and working-tree marker.
%
% Return unavailable if Git cannot resolve HEAD; append +uncommitted when
% the working-tree query succeeds and reports changes.

[status, revision] = system('git rev-parse --short HEAD');
if status == 0
    code_version = string(strtrim(revision));
    [dirty_status, changes] = system('git status --porcelain --untracked-files=all');
    if dirty_status == 0 && strlength(strtrim(string(changes))) > 0
        code_version = code_version + "+uncommitted";
    end
else
    code_version = "unavailable";
end
end
