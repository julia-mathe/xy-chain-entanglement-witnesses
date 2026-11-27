function [c, mz] = spin_correlations_symbol_tdl(Mhat, R, k)
% SPIN_CORRELATIONS_SYMBOL_TDL  Connected two-spin correlations of a translation-invariant Gaussian state.

if nargin < 3 || isempty(k)
    k = momentum_grid_tdl();
end
Mw = majorana_window_tdl(Mhat, R + 1, k);
mz = -Mw(1, 2);
c = zeros(R + 1, 3, 3);
c(1, :, :) = reshape(diag([1, 1, 1 - mz^2]), 1, 3, 3);
for r = 1:R
    for a = 1:3
        for b = 1:3
            if xor(a == 3, b == 3)
                continue
            end
            v = real(pauli_string_expectation(Mw, [1, r + 1], [a, b]));
            if a == 3 && b == 3
                v = v - mz^2;
            end
            c(r + 1, a, b) = v;
        end
    end
end
end
