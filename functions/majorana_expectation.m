function value = majorana_expectation(M, idx)
% MAJORANA_EXPECTATION  <a_m1 ... a_m2n> = (-i)^n Pf(M_m) in a Gaussian state (App. G).

n = numel(idx);
if mod(n, 2) == 1
    value = 0;
    return
end
K = zeros(n);
for p = 1:n
    for q = p+1:n
        K(p, q) = (idx(p) == idx(q)) - 1i * M(idx(p), idx(q));
        K(q, p) = -K(p, q);
    end
end
value = pfaffian_parlett_reid(K);
end
