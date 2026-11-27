function X = partial_transpose(rho, A, N)
% PARTIAL_TRANSPOSE  Partial transpose of an N-spin matrix on the spins in A.

R = reshape(rho, 2*ones(1, 2*N));
order = 1:2*N;
for k = A(:)'
    row_dimension    = N - k + 1;
    column_dimension = 2*N - k + 1;
    order([row_dimension, column_dimension]) = order([column_dimension, row_dimension]);
end
X = reshape(permute(R, order), 2^N, 2^N);
end
