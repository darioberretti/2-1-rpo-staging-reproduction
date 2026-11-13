function save_result(result, filename)
%SAVE_RESULT Save one reusable result using the documented top-level variable.
%
% PROTOTYPE:
%   rpo.save_result(result, filename)
%
% DESCRIPTION:
%   Create the parent directory if needed and write result in MAT v7.3
%   format. An existing file is replaced; overwrite policy belongs to the
%   caller, normally rpo.prepare_run.
%
% INPUT:
%   result    Study result to store as the top-level variable result.
%   filename  Destination MAT file path.
%
% OUTPUT:
%   None. Writes the result MAT file.

folder = fileparts(filename);
if ~isfolder(folder)
    mkdir(folder);
end
save(filename, 'result', '-v7.3');
end
