function [omega, s, qmin, best, S, q] = cmc_theorem1_ring(Gamma)
% CMC_THEOREM1_RING  omega_CMC of a ring from Theorem 1; minimizing momenta folded to [0, pi].

[q, S] = structure_factors_ring(Gamma);
[s, k] = min(S, [], 2);
s = s';
qmin = q(k);
qmin = min(qmin, 2*pi - qmin);
[omega, best] = theorem1_omega(s);
end
