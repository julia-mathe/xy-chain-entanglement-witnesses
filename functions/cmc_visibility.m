function p_max = cmc_visibility(m, C, N, periodic)
% CMC_VISIBILITY  Largest noise factor p for which the CMC detects the rescaled correlations.

low = 0;
high = 1.5;
if omega_of_shrunk_data(high, m, C, N, periodic) >= 1
    p_max = Inf;
    return
end
for halving = 1:40
    middle = (low + high) / 2;
    if omega_of_shrunk_data(middle, m, C, N, periodic) >= 1
        low = middle;
    else
        high = middle;
    end
end
p_max = (low + high) / 2;
end

function omega = omega_of_shrunk_data(p, m, C, N, periodic)
mp = p * m;
Cp = p * C;
for j = 1:N
    block = 3*j-2 : 3*j;
    Cp(block, block) = eye(3);
end
Gamma = Cp - mp * mp';
if periodic
    omega = cmc_theorem1_ring(Gamma);
else
    omega = cmc_sdp(Gamma);
end
end
