function e = energy_product_state_min(gamma, h)
% ENERGY_PRODUCT_STATE_MIN  Minimum energy per site over product states (App. A).

a = (1 + gamma) / 2;
if h < 2*a
    e = -a - h^2 / (4*a);
else
    e = -h;
end
end
