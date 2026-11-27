function rho = reduced_density_matrix_gaussian(M, sites)
% REDUCED_DENSITY_MATRIX_GAUSSIAN  Reduced state of a few spins of a Gaussian state.

L = numel(sites);
[sx, sy, sz] = pauli_matrices();
paulis = {speye(2), sx, sy, sz};
rho = zeros(2^L);
for s = 0:4^L - 1
    labels = zeros(1, L);
    rest = s;
    for p = L:-1:1
        labels(p) = mod(rest, 4);
        rest = floor(rest / 4);
    end
    on = labels > 0;
    if any(on)
        value = pauli_string_expectation(M, sites(on), labels(on));
    else
        value = 1;
    end
    op = 1;
    for p = 1:L
        op = kron(op, full(paulis{labels(p) + 1}));
    end
    rho = rho + value * op;
end
rho = rho / 2^L;
rho = (rho + rho') / 2;
end
