function key = monomial_key(monomial)
% MONOMIAL_KEY  Text key of a monomial (dictionary lookup).

key = ['m', sprintf('_%d', sort(monomial))];
end
