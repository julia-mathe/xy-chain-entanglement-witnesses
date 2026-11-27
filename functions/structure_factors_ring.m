function [q, S] = structure_factors_ring(Gamma)
% STRUCTURE_FACTORS_RING  Structure factors of a ring, Eq. (16).

N = size(Gamma, 1) / 3;
q = 2*pi*(0:N-1) / N;
S = zeros(3, N);
sites = (1:N)';

for alpha = 1:3
    idx = 3*(sites - 1) + alpha;
    G = Gamma(idx, idx);
    for k = 1:N
        v = exp(1i * q(k) * sites) / sqrt(N);
        S(alpha, k) = real(v' * G * v);
    end
end
end
