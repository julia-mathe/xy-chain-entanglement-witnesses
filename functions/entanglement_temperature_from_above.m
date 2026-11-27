function T_ent = entanglement_temperature_from_above(omega_of_T, T_high, T_low, number_of_T)
% ENTANGLEMENT_TEMPERATURE_FROM_ABOVE  T_ent found by lowering T step by step.

temperatures = logspace(log10(T_high), log10(T_low), number_of_T);
all_undetected = true;
previous_T = NaN;
previous_value = NaN;

for n = 1:number_of_T
    T = temperatures(n);
    value = omega_of_T(T);
    if value < 1
        all_undetected = false;
    end
    if n > 1 && previous_value >= 1 && value < 1
        T_ent = fzero(@(t) omega_of_T(t) - 1, [T, previous_T], optimset('TolX', 1e-4));
        return
    end
    previous_T = T;
    previous_value = value;
end

if all_undetected
    T_ent = 0;
else
    T_ent = NaN;
end
end
