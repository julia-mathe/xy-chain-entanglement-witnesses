function [Gamma, m, C] = covariance_matrix(rho, ops)
% COVARIANCE_MATRIX  Covariance matrix of all local Pauli operators, Eq. (7).

N = size(ops, 1);
M = 3 * N;

A = cell(M, 1);
for j = 1:N
    for alpha = 1:3
        A{3*(j-1) + alpha} = ops{j, alpha};
    end
end

m = zeros(M, 1);
for p = 1:M
    m(p) = expectation_value(rho, A{p});
end

C = zeros(M, M);
for p = 1:M
    for q = p:M
        C(p, q) = expectation_value(rho, A{p}*A{q} + A{q}*A{p}) / 2;
        C(q, p) = C(p, q);
    end
end

Gamma = C - m * m';
end
