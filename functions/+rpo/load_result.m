function result = load_result(filename)
%LOAD_RESULT Load the documented top-level result variable from a MAT file.
%
% PROTOTYPE:
%   result = rpo.load_result(filename)
%
% DESCRIPTION:
%   Require an existing MAT file containing a variable named result. This
%   checks the container, leaving study-specific fields to its consumers.
%
% INPUT:
%   filename  Path to a reference or generated result MAT file.
%
% OUTPUT:
%   result    Stored result variable; run01 consumers expect a study struct.

if ~isfile(filename)
    error('rpo:ResultMissing', [ ...
        'Required result file is not installed or has not been generated:\n%s\n\n', ...
        'For reference results, extract accepted_results.zip into the repository root ', ...
        '(see README.md). For local results, run the experiment and save them.'], filename);
end
loaded = load(filename, 'result');
if ~isfield(loaded, 'result')
    error('rpo:ResultSchema', ...
        'The MAT file does not contain the required top-level variable "result":\n%s', ...
        filename);
end
result = loaded.result;
end
