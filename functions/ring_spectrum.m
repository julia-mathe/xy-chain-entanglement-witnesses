function ring = ring_spectrum(N, gamma, h)
% RING_SPECTRUM  Spectrum and correlations of a ring, for evaluating many temperatures.

ops = local_pauli_operators(N);
H = xy_hamiltonian(N, gamma, h, true);
[V, E] = eig(full(H));

ring.N = N;
ring.gamma = gamma;
ring.h = h;
ring.E = diag(E);
ring.rmax = floor(N/2);

number_of_states = numel(ring.E);
ring.corr = zeros(number_of_states, ring.rmax + 1, 3);
for alpha = 1:3
    ring.corr(:, 1, alpha) = 1;
    for r = 1:ring.rmax
        X = ops{1, alpha} * ops{1 + r, alpha};
        ring.corr(:, r + 1, alpha) = real(sum(conj(V) .* (X * V), 1)).';
    end
end
ring.mz = real(sum(conj(V) .* (ops{1, 3} * V), 1)).';
end
