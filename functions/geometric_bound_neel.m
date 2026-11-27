function [bound, mu_best] = geometric_bound_neel(omega, N)
% GEOMETRIC_BOUND_NEEL  Lower bound on E_G from Eq. (36) with the inner maximum at Neel-type states.

[ells, p] = neel_total_spin_distribution(N);
w = 2*ells.*(ells + 1)/N - 1;

g = @(mu) largest_secular_root(p, -mu * w);
minus_bound = @(mu) -(1 + mu*(1 - omega) - g(mu));

[mu_best, value] = fminbnd(minus_bound, 1e-6, 20, optimset('TolX', 1e-10));
bound = -value;
end
