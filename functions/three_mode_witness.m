function w = three_mode_witness(Gamma, q)
% THREE_MODE_WITNESS  omega_q of Eq. (19) for momenta q in {0, pi}.

N = size(Gamma, 1) / 3;
w = 0;
for alpha = 1:3
    u = zeros(3*N, 1);
    for j = 1:N
        u(3*(j-1) + alpha) = cos(q(alpha) * (j-1)) / sqrt(N);
    end
    w = w + u' * Gamma * u;
end
w = w / 2;
end
