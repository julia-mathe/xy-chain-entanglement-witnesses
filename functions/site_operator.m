function O = site_operator(op, j, N)
% SITE_OPERATOR  Single-spin operator op acting on spin j of N spins.

left  = speye(2^(j-1));
right = speye(2^(N-j));
O = kron(kron(left, op), right);
end
