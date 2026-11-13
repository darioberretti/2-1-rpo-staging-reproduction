function [c, ceq, Jc, Jceq] = constr_CR3BP_reachability(yy, problem)
%CONSTR_CR3BP_REACHABILITY Dense fixed-phase reachability constraints.
%
% PROTOTYPE:
%   [c, ceq, Jc, Jceq] = constr_CR3BP_reachability(yy, problem)
%
% DESCRIPTION:
%   Set the objective-mode/phase flags and delegate constraint evaluation
%   to constr_CR3BP_phase_core. Keep the shared dense Jacobian storage.
%
% INPUT:
%   yy       Decision vector from RPO2DS_transfer setup.
%   problem  Transfer problem_data struct or supported legacy constants; see
%              constr_CR3BP_phase_core for layout, units, and active policies.
%
% OUTPUT:
%   c, ceq    Inequality (<= 0) and equality (= 0) residual columns.
%   Jc, Jceq  Optional [numel(yy) x number of constraints] Jacobians.

[c, ceq, Jc, Jceq] = constr_CR3BP_phase_core(yy, problem);
end
