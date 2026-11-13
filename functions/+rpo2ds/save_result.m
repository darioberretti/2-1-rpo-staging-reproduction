function save_result(result, filename, overwrite_allowed)
%SAVE_RESULT Save one compact Run 02 result with explicit overwrite policy.
%
% PROTOTYPE:
%   rpo2ds.save_result(result, filename, overwrite_allowed)
%
% DESCRIPTION:
%   Create the parent directory and write result in MAT v7.3 format. Refuse
%   to replace an existing file unless overwrite_allowed is true.
%
% INPUT:
%   result             Compact study result struct.
%   filename           Destination MAT file path.
%   overwrite_allowed  Optional scalar logical; defaults to false.
%
% OUTPUT:
%   None. Writes the top-level variable result.

if nargin < 3
    overwrite_allowed = false;
end
if isfile(filename) && ~overwrite_allowed
    error('rpo2ds:ResultExists', ...
        'Refusing to replace an existing result without overwrite permission:\n%s', filename);
end
parent = fileparts(filename);
if ~isfolder(parent)
    mkdir(parent);
end
save(filename, 'result', '-v7.3');
end
