function [C, ells, D, report] = ppt_weight_constraints(N)
% PPT_WEIGHT_CONSTRAINTS  Linear PPT conditions on the total-spin weights (App. H).

[ells, P, D] = total_spin_projectors(N);
L = numel(ells);
C = zeros(0, L);
report = zeros(N/2, 2);
rng(0);

for k = 1:N/2
    A = 1:k;
    X = cell(1, L);
    for n = 1:L
        X{n} = partial_transpose(P{n}, A, N);
    end
    c = randn(1, L);
    mixture = zeros(2^N);
    for n = 1:L
        mixture = mixture + c(n) * X{n};
    end
    [U, ~] = eig((mixture + mixture') / 2);
    values = zeros(L, 2^N);
    largest_offdiagonal = 0;
    for n = 1:L
        Y = U' * X{n} * U;
        values(n, :) = diag(Y)';
        largest_offdiagonal = max(largest_offdiagonal, max(max(abs(Y - diag(diag(Y))))));
    end
    rows = values' ./ D;
    rows = unique(round(rows * 1e10) / 1e10, 'rows');
    C = [C; rows];
    report(k, :) = [largest_offdiagonal, size(rows, 1)];
end
end
