function X = partial_transpose_by_loops(rho, A, N)
% PARTIAL_TRANSPOSE_BY_LOOPS  Slow reference implementation of the partial transpose.

dimension = 2^N;
X = zeros(dimension);
powers = 2.^(N-1:-1:0);

for i = 0:dimension-1
    bits_i = bitget(i, N:-1:1);
    for j = 0:dimension-1
        bits_j = bitget(j, N:-1:1);
        new_i = bits_i;
        new_j = bits_j;
        new_i(A) = bits_j(A);
        new_j(A) = bits_i(A);
        X(sum(new_i .* powers) + 1, sum(new_j .* powers) + 1) = rho(i + 1, j + 1);
    end
end
end
