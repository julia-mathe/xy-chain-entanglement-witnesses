function TE = energy_witness_temperature(gamma, h)
% ENERGY_WITNESS_TEMPERATURE  Highest temperature at which the energy witness detects entanglement.

difference = @(T) energy_per_site_tdl(gamma, h, T) - energy_product_state_min(gamma, h);

if difference(1e-3) >= 0
    TE = 0;
    return
end
TE = fzero(difference, [1e-3, 10], optimset('TolX', 1e-6));
end
