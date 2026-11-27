function v = moment_hierarchy_value(rho, N, k, solver)
% MOMENT_HIERARCHY_VALUE  Optimal value v_SDP of the moment hierarchy at order k (Eq. (15), App. C).
%   All 4^N Pauli correlations of rho are the data; v < 0 certifies entanglement.
if nargin < 4, solver = ''; end
require_cvx();
n = 3*N;
MI = multiindices(n, 2*k);
dict = containers.Map('KeyType', 'char', 'ValueType', 'double');
for j = 1:size(MI, 1)
    dict(sprintf('%d_', MI(j, :))) = j;
end
Bk  = multiindices(n, k);
Bkm = multiindices(n, k - 1);
[known, values] = pauli_moments_all(rho, N);
known_idx = zeros(size(values));
for j = 1:numel(values)
    known_idx(j) = dict(sprintf('%d_', known(j, :)));
end
free = true(size(MI, 1), 1);
free(known_idx) = false;
I0 = zeros(size(Bk, 1));
for a = 1:size(Bk, 1)
    for b = 1:size(Bk, 1)
        I0(a, b) = dict(sprintf('%d_', Bk(a, :) + Bk(b, :)));
    end
end
nl = size(Bkm, 1);
Iloc = zeros(nl, nl, N, 4);
for q = 1:N
    for a = 1:nl
        for b = 1:nl
            base = Bkm(a, :) + Bkm(b, :);
            Iloc(a, b, q, 1) = dict(sprintf('%d_', base));
            for alpha = 1:3
                e = zeros(1, n);
                e(3*(q-1) + alpha) = 2;
                Iloc(a, b, q, alpha + 1) = dict(sprintf('%d_', base + e));
            end
        end
    end
end
cvx_begin sdp quiet
    if ~isempty(solver), cvx_solver(solver); end
    variable Q0(size(Bk, 1), size(Bk, 1)) symmetric
    variable Qloc(nl, nl, N) symmetric
    variable c(numel(values))
    expression coeff(size(MI, 1))
    for a = 1:size(Bk, 1)
        for b = 1:size(Bk, 1)
            coeff(I0(a, b)) = coeff(I0(a, b)) + Q0(a, b);
        end
    end
    for q = 1:N
        for a = 1:nl
            for b = 1:nl
                coeff(Iloc(a, b, q, 1)) = coeff(Iloc(a, b, q, 1)) + Qloc(a, b, q);
                for alpha = 2:4
                    coeff(Iloc(a, b, q, alpha)) = coeff(Iloc(a, b, q, alpha)) - Qloc(a, b, q);
                end
            end
        end
    end
    normalization = trace(Q0);
    for q = 1:N
        normalization = normalization + trace(Qloc(:, :, q));
    end
    minimize( values.' * c )
    subject to
        Q0 >= 0;
        for q = 1:N
            Qloc(:, :, q) >= 0;
        end
        normalization <= 1;
        coeff(known_idx) == c;
        coeff(free) == 0;
cvx_end
v = cvx_optval;
end

function MI = multiindices(n, dmax)
% All exponent vectors of n variables with total degree <= dmax.
if n == 1
    MI = (0:dmax)';
    return
end
MI = zeros(0, n);
for v = 0:dmax
    rest = multiindices(n - 1, dmax - v);
    MI = [MI; v * ones(size(rest, 1), 1), rest]; %#ok<AGROW>
end
end
