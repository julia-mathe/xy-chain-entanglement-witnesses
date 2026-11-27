function d = leading_principal_minors(M, stop)
% LEADING_PRINCIPAL_MINORS  All leading principal minors (Toeplitz determinants, Eqs. (43)-(44)).

if nargin < 2
    stop = 1e-15;
end

n = size(M, 1);
d = zeros(1, n);
A = M;
product_of_pivots = 1;

for k = 1:n
    pivot = A(k, k);
    product_of_pivots = product_of_pivots * pivot;
    d(k) = product_of_pivots;
    if abs(product_of_pivots) < stop || abs(pivot) < 1e-300
        break
    end
    if k < n
        factors = A(k+1:n, k) / pivot;
        A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - factors * A(k, k+1:n);
    end
end
end
