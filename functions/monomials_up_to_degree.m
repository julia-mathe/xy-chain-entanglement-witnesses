function monomials = monomials_up_to_degree(number_of_variables, max_degree)
% MONOMIALS_UP_TO_DEGREE  All monomials of degree <= max_degree.

n = number_of_variables;
monomials = {[]};
for d = 1:max_degree
    increasing = nchoosek(1:n + d - 1, d);
    nondecreasing = increasing - (0:d-1);
    for row = 1:size(nondecreasing, 1)
        monomials{end + 1} = nondecreasing(row, :);
    end
end
end
