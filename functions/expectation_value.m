function value = expectation_value(rho, X)
% EXPECTATION_VALUE  Real part of Tr(rho X).

value = real(sum(sum(rho.' .* X)));
end
