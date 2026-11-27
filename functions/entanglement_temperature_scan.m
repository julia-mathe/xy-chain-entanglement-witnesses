function [T_ent, temperatures, values] = entanglement_temperature_scan(omega_of_T, T_low, T_high, number_of_T)
% ENTANGLEMENT_TEMPERATURE_SCAN  T_ent = highest temperature with omega_CMC < 1.

temperatures = logspace(log10(T_low), log10(T_high), number_of_T);
values = zeros(size(temperatures));
for n = 1:number_of_T
    values(n) = omega_of_T(temperatures(n));
end

detected = values < 1;
if ~any(detected)
    T_ent = 0;
    return
end
last = find(detected, 1, 'last');
if last == number_of_T
    T_ent = Inf;
    return
end
T_ent = fzero(@(T) omega_of_T(T) - 1, [temperatures(last), temperatures(last + 1)], ...
              optimset('TolX', 1e-4));
end
