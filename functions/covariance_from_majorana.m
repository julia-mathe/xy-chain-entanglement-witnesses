function [Gamma, mz] = covariance_from_majorana(M, sites)
% COVARIANCE_FROM_MAJORANA  Spin covariance matrix of selected spins of a fermionic Gaussian state.

L = numel(sites);
mz = zeros(L, 1);
for p = 1:L
    j = sites(p);
    mz(p) = -M(2*j-1, 2*j);
end

Gamma = zeros(3*L);
for p = 1:L
    Gamma(3*p-2, 3*p-2) = 1;
    Gamma(3*p-1, 3*p-1) = 1;
    Gamma(3*p,   3*p)   = 1 - mz(p)^2;
    for q = p+1:L
        i = sites(p);  j = sites(q);
        for a = 1:3
            for b = 1:3
                if i < j
                    v = real(pauli_string_expectation(M, [i j], [a b]));
                else
                    v = real(pauli_string_expectation(M, [j i], [b a]));
                end
                if a == 3 && b == 3
                    v = v - mz(p) * mz(q);
                end
                Gamma(3*(p-1)+a, 3*(q-1)+b) = v;
                Gamma(3*(q-1)+b, 3*(p-1)+a) = v;
            end
        end
    end
end
end
