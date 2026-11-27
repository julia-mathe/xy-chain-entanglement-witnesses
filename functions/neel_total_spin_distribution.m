function [ells, p] = neel_total_spin_distribution(N)
% NEEL_TOTAL_SPIN_DISTRIBUTION  Total-spin weights of a Neel state.

n = N / 2;
ells = 0:n;
logp = log(2*ells + 1) + 2*gammaln(n + 1) - gammaln(n + ells + 2) - gammaln(n - ells + 1);
p = exp(logp);
end
