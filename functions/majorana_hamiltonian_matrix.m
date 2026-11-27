function A = majorana_hamiltonian_matrix(N, gamma, h, J)
% MAJORANA_HAMILTONIAN_MATRIX  Matrix A of H = (i/4) a^T A a for an open chain (App. G).

if nargin < 4
    J = 1;
end

A = zeros(2*N);
for j = 1:N
    A(2*j-1, 2*j) = 2*h;
end
for j = 1:N-1
    A(2*j, 2*j+1)   =  J*(1 + gamma);
    A(2*j-1, 2*j+2) = -J*(1 - gamma);
end
A = A - A.';
end
