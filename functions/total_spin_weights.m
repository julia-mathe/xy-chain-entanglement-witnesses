function p = total_spin_weights(psi, V, ell_of_state, ells)
% TOTAL_SPIN_WEIGHTS  Weights of a pure state in the total-spin sectors.

amplitudes = abs(V' * psi).^2;
p = zeros(size(ells));
for n = 1:numel(ells)
    p(n) = sum(amplitudes(ell_of_state == ells(n)));
end
end
