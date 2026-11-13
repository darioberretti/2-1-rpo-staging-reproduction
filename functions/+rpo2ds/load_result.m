function result = load_result(filename)
%LOAD_RESULT Load one Run 02 result with a useful missing-file error.
%
% PROTOTYPE:
%   result = rpo2ds.load_result(filename)
%
% DESCRIPTION:
%   Require an existing MAT file containing a top-level result variable.
%   Study-specific fields remain the responsibility of the consuming function.
%
% INPUT:
%   filename  Reference or generated result MAT file path.
%
% OUTPUT:
%   result    Stored variable containing the compact study result.

if ~isfile(filename)
    error('rpo2ds:ResultMissing', [ ...
        'Required Run 02 result is missing:\n%s\n\n', ...
        'For reference results, extract accepted_results.zip into the repository root ', ...
        '(see README.md). For local results, run Run 02 with workflow = "run_experiment" ', ...
        'and save_computed_results = true.'], filename);
end
loaded = load(filename, 'result');
if ~isfield(loaded, 'result') || ~isstruct(loaded.result)
    error('rpo2ds:ResultSchema', ...
        'The MAT file must contain one top-level struct named "result":\n%s', filename);
end
result = loaded.result;
required = {'data', 'metrics', 'configuration', 'provenance', 'status'};
missing = required(~isfield(result, required));
if ~isempty(missing)
    error('rpo2ds:ResultSchema', ...
        'Run 02 result is missing required fields (%s):\n%s', ...
        strjoin(missing, ', '), filename);
end
end
