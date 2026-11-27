function omega = omega_structure_factors_tdl(c, nq)
% OMEGA_STRUCTURE_FACTORS_TDL  Right-hand side of Eq. (25) for a translation-invariant state.

if nargin < 2
    nq = [];
end
[~, ~, omega] = structure_factor_minima_tdl(c, nq);
end
