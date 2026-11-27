function [sx, sy, sz] = pauli_matrices()
% PAULI_MATRICES  The Pauli matrices (sparse).

sx = sparse([0 1; 1 0]);
sy = sparse([0 -1i; 1i 0]);
sz = sparse([1 0; 0 -1]);
end
