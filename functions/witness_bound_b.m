function b = witness_bound_b(Z)
% WITNESS_BOUND_B  Separable bound b(Z), Eq. (10).

N = size(Z, 1) / 3;
b = 0;
for j = 1:N
    block = 3*j-2 : 3*j;
    Zjj = (Z(block, block) + Z(block, block)') / 2;
    b = b + trace(Zjj) - max(eig(Zjj));
end
end
