function [Ahat, k] = majorana_symbol_tdl(gamma, h, k, J)
% MAJORANA_SYMBOL_TDL  2 x 2 momentum-space blocks of the Majorana Hamiltonian.

if nargin < 3 || isempty(k)
    k = momentum_grid_tdl();
end
if nargin < 4
    J = 1;
end
A0  = [0, 2*h; -2*h, 0];
A1  = [0, -J*(1 - gamma); J*(1 + gamma), 0];
Am1 = -A1.';
k = k(:).';
phase = reshape(exp(-1i * k), 1, 1, []);
Ahat = A0 + A1 .* phase + Am1 .* conj(phase);
end
