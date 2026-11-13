function tf = RPO2DS_has_patch_solution(FullSolution)
%RPO2DS_HAS_PATCH_SOLUTION Test whether the patch output contains usable selected rows.
%
% PROTOTYPE:
%   has_solution = RPO2DS_has_patch_solution(solutionTable)
%
% DESCRIPTION:
%   Require the expected table fields and finite perigee cost used downstream.
%
% INPUT:
%   solutionTable  Patch-selection table or an empty/unsupported value.
%
% OUTPUT:
%   has_solution   Scalar logical; see the implemented required-field checks.

tf = istable(FullSolution) && height(FullSolution) > 0 && ...
    ismember('deltavpMag', FullSolution.Properties.VariableNames) && ...
    ~isempty(FullSolution.deltavpMag) && isfinite(FullSolution.deltavpMag(1));
end
