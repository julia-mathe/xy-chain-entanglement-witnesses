function eta = largest_secular_root(p, a)
% LARGEST_SECULAR_ROOT  Largest eigenvalue of lambda*W' + |phi><phi| from the total-spin weights of phi.

used = p(:) > 1e-15;
v  = sqrt(p(used));  v  = v(:);
aa = a(used);        aa = aa(:);
eta = max(eig(diag(aa) + v * v'));
eta = max(eta, max(a));
end
