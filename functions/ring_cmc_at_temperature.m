function [omega, s, qmin, S, q, best] = ring_cmc_at_temperature(ring, T)
% RING_CMC_AT_TEMPERATURE  omega_CMC of a ring at temperature T from Theorem 1.

N = ring.N;
E = ring.E;

if T > 0
    w = exp(-(E - min(E)) / T);
else
    w = double((1:numel(E))' == 1);
end
w = w / sum(w);

mz = w' * ring.mz;
q = 2*pi*(0:N-1) / N;
r = 0:N-1;
S = zeros(3, N);

for alpha = 1:3
    Cr = w' * ring.corr(:, :, alpha);
    if alpha == 3
        Cr = Cr - mz^2;
    end
    Cfull = Cr(min(r, N - r) + 1);
    for k = 1:N
        S(alpha, k) = real(sum(Cfull .* exp(1i * q(k) * r)));
    end
end

[s, kmin] = min(S, [], 2);
s = s';
qmin = q(kmin);
[omega, best] = theorem1_omega(s);
end
