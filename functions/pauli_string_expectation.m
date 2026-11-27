function value = pauli_string_expectation(M, sites, components)
% PAULI_STRING_EXPECTATION  <sigma_j1^a1 ... sigma_jn^an> in a Gaussian state (Jordan-Wigner, App. G).

idx = [];
coeff = 1;
for n = 1:numel(sites)
    j = sites(n);
    if components(n) == 3
        idx = [idx, 2*j-1, 2*j];
        coeff = coeff * (-1i);
    else
        idx = [idx, 1:2*(j-1)];
        coeff = coeff * (-1i)^(j-1);
        if components(n) == 1
            idx = [idx, 2*j-1];
        else
            idx = [idx, 2*j];
        end
    end
end

inversions = sum(sum(triu(idx(:) > idx(:).', 1)));
sign_factor = (-1)^inversions;

[u, ~, k] = unique(idx);
multiplicity = accumarray(k(:), 1);
remaining = u(mod(multiplicity, 2) == 1);

if mod(numel(remaining), 2) == 1
    value = 0;
elseif isempty(remaining)
    value = coeff * sign_factor;
else
    value = coeff * sign_factor * majorana_expectation(M, remaining);
end
end
