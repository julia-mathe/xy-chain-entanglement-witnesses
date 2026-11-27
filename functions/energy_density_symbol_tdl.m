function e = energy_density_symbol_tdl(Mhat, gamma, h, k)
% ENERGY_DENSITY_SYMBOL_TDL  Energy per site of a translation-invariant Gaussian state.

if nargin < 4 || isempty(k)
    k = momentum_grid_tdl();
end
Ahat = majorana_symbol_tdl(gamma, h, k);
total = 0;
for n = 1:numel(k)
    total = total + trace(Ahat(:, :, n) * Mhat(:, :, n));
end
e = real(-total / (4 * numel(k)));
end
