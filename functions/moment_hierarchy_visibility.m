function p_max = moment_hierarchy_visibility(m, C, N, k)
% MOMENT_HIERARCHY_VISIBILITY  Moment hierarchy with one- and two-spin data: largest detected noise factor.

require_cvx();
number_of_variables = 3 * N;
basis = monomials_up_to_degree(number_of_variables, k);
all_monomials = monomials_up_to_degree(number_of_variables, 2*k);
B = numel(basis);
number_of_moments = numel(all_monomials);

position = containers.Map();
for n = 1:number_of_moments
    position(monomial_key(all_monomials{n})) = n;
end
position_of = @(mono) position(monomial_key(mono));

rows = zeros(B*B, 1);
columns = zeros(B*B, 1);
count = 0;
for a = 1:B
    for b = 1:B
        count = count + 1;
        rows(count) = a + (b - 1)*B;
        columns(count) = position_of([basis{a}, basis{b}]);
    end
end
Selection = sparse(rows, columns, 1, B*B, number_of_moments);

index_one = position_of([]);
index1 = zeros(number_of_variables, 1);
value1 = zeros(number_of_variables, 1);
for v = 1:number_of_variables
    index1(v) = position_of(v);
    value1(v) = m(v);
end
index2 = [];
value2 = [];
for i = 1:N
    for j = i+1:N
        for a = 1:3
            for b = 1:3
                vi = 3*(i-1) + a;
                vj = 3*(j-1) + b;
                index2(end+1, 1) = position_of([vi, vj]);
                value2(end+1, 1) = C(vi, vj);
            end
        end
    end
end

low_monomials = monomials_up_to_degree(number_of_variables, 2*k - 2);
index0 = [];  indexx = [];  indexy = [];  indexz = [];
for n = 1:numel(low_monomials)
    mono = low_monomials{n};
    for i = 1:N
        ix = 3*(i-1) + 1;   iy = 3*(i-1) + 2;   iz = 3*(i-1) + 3;
        index0(end+1, 1) = position_of(mono);
        indexx(end+1, 1) = position_of([mono, ix, ix]);
        indexy(end+1, 1) = position_of([mono, iy, iy]);
        indexz(end+1, 1) = position_of([mono, iz, iz]);
    end
end

cvx_begin sdp quiet
    variable y(number_of_moments)
    variable p
    maximize( p )
    subject to
        reshape(Selection * y, B, B) >= 0;
        y(index_one) == 1;
        y(index1) == p * value1;
        y(index2) == p * value2;
        y(indexx) + y(indexy) + y(indexz) == y(index0);
        p <= 5;
cvx_end

p_max = p;
end
