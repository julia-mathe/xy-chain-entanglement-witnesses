function ops = local_pauli_operators(N)
% LOCAL_PAULI_OPERATORS  All sigma_j^alpha of N spins (spin 1 is the leftmost tensor factor).

[sx, sy, sz] = pauli_matrices();
single_spin = {sx, sy, sz};

ops = cell(N, 3);
for j = 1:N
    for alpha = 1:3
        ops{j, alpha} = site_operator(single_spin{alpha}, j, N);
    end
end
end
