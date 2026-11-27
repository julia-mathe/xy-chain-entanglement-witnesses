function H = xy_hamiltonian(N, gamma, h, periodic, J)
% XY_HAMILTONIAN  XY chain in a transverse field, Eq. (1).

if nargin < 5
    J = 1;
end

ops = local_pauli_operators(N);

if periodic
    number_of_bonds = N;
else
    number_of_bonds = N - 1;
end

H = sparse(2^N, 2^N);
for j = 1:number_of_bonds
    k = mod(j, N) + 1;
    H = H - J * ( (1+gamma)/2 * ops{j,1} * ops{k,1} ...
                + (1-gamma)/2 * ops{j,2} * ops{k,2} );
end
for j = 1:N
    H = H - h * ops{j,3};
end

H = real(H);
end
