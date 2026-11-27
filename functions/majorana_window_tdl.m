function Mw = majorana_window_tdl(Mhat, L, k)
% MAJORANA_WINDOW_TDL  Majorana correlation matrix of L consecutive spins of the infinite chain.

if nargin < 3 || isempty(k)
    k = momentum_grid_tdl();
end
nk = numel(k);
dmax = L - 1;
Md = zeros(2, 2, 2*dmax + 1);
for d = -dmax:dmax
    ph = reshape(exp(1i * k * d), 1, 1, nk);
    Md(:, :, d + dmax + 1) = real(sum(Mhat .* ph, 3)) / nk;
end
Mw = zeros(2*L);
for j = 1:L
    for l = 1:L
        Mw(2*j-1:2*j, 2*l-1:2*l) = Md(:, :, (l - j) + dmax + 1);
    end
end
Mw = (Mw - Mw.') / 2;
end
