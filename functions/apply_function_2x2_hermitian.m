function [F, P_plus, P_minus, degenerate] = apply_function_2x2_hermitian(X, f)
% APPLY_FUNCTION_2X2_HERMITIAN  f(X) for a stack of Hermitian 2 x 2 matrices X(:,:,k).

n = size(X, 3);
a  = real(reshape(X(1, 1, :), [], 1));
d  = real(reshape(X(2, 2, :), [], 1));
b  = reshape(X(1, 2, :), [], 1);
x0 = (a + d) / 2;
xz = (a - d) / 2;
r  = sqrt(xz.^2 + abs(b).^2);
degenerate = r < 1e-14;
r(degenerate) = 0;
rs = r;
rs(degenerate) = 1;

Yn = zeros(2, 2, n);
Yn(1, 1, :) = xz ./ rs;
Yn(2, 2, :) = -xz ./ rs;
Yn(1, 2, :) = b ./ rs;
Yn(2, 1, :) = conj(b) ./ rs;
Yn(:, :, degenerate) = 0;

I2 = repmat(eye(2), [1, 1, n]);
P_plus  = (I2 + Yn) / 2;
P_minus = (I2 - Yn) / 2;
fp = reshape(f(x0 + r), 1, 1, n);
fm = reshape(f(x0 - r), 1, 1, n);
F = P_plus .* fp + P_minus .* fm;
end
