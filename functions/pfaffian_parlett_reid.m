function pf = pfaffian_parlett_reid(K)
% PFAFFIAN_PARLETT_REID  Pfaffian of an antisymmetric matrix.

n = size(K, 1);
if mod(n, 2) == 1
    pf = 0;
    return
end
A = K;
pf = 1;
for k = 1:2:n-1
    [~, p] = max(abs(A(k+1:n, k)));
    kp = k + p;
    if kp ~= k + 1
        A([k+1, kp], :) = A([kp, k+1], :);
        A(:, [k+1, kp]) = A(:, [kp, k+1]);
        pf = -pf;
    end
    if A(k+1, k) == 0
        pf = 0;
        return
    end
    pf = pf * A(k, k+1);
    if k + 2 <= n
        tau = A(k, k+2:n) / A(k, k+1);
        u = A(k+2:n, k+1);
        A(k+2:n, k+2:n) = A(k+2:n, k+2:n) + tau.' * u.' - u * tau;
    end
end
end
