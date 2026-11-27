function psi = product_state(theta, phi)
% PRODUCT_STATE  Pure product state of N spins from Bloch angles.

psi = 1;
for j = 1:numel(theta)
    single_spin = [cos(theta(j)/2); exp(1i*phi(j)) * sin(theta(j)/2)];
    psi = kron(psi, single_spin);
end
end
