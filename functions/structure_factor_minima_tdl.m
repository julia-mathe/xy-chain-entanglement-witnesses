function [s, qmin, omega, q, S] = structure_factor_minima_tdl(c, nq)
% STRUCTURE_FACTOR_MINIMA_TDL  Minimal structure factors, Eq. (17), and Eq. (25) for the infinite chain.

if nargin < 2 || isempty(nq)
    nq = 2001;
end
q = linspace(0, pi, nq);
Ghat = covariance_fourier_tdl(c, q);
S = zeros(3, nq);
for a = 1:3
    S(a, :) = real(squeeze(Ghat(a, a, :))).';
end
[s, kmin] = min(S, [], 2);
s = s.';
qmin = q(kmin);
omega = theorem1_omega(s);
end
