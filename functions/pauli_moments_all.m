function [alphas, values] = pauli_moments_all(rho, N)
% PAULI_MOMENTS_ALL  All 4^N Pauli correlations <sigma^mu_1 ... sigma^mu_N> and their exponent vectors.
P = {eye(2), [0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1]};
alphas = zeros(4^N, 3*N);
values = zeros(4^N, 1);
for s = 0:4^N - 1
    mu = zeros(1, N);
    t = s;
    for q = N:-1:1
        mu(q) = mod(t, 4);
        t = floor(t / 4);
    end
    op = 1;
    for q = 1:N
        op = kron(op, P{mu(q) + 1});
        if mu(q) > 0
            alphas(s + 1, 3*(q-1) + mu(q)) = 1;
        end
    end
    values(s + 1) = real(trace(rho * op));
end
end
