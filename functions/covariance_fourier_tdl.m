function Ghat = covariance_fourier_tdl(c, q)
% COVARIANCE_FOURIER_TDL  3 x 3 Fourier blocks of the covariance matrix of a translation-invariant state.

R = size(c, 1) - 1;
Ghat = zeros(3, 3, numel(q));
for n = 1:numel(q)
    G = squeeze(c(1, :, :));
    for d = 1:R
        cd = squeeze(c(d + 1, :, :));
        G = G + exp(-1i * q(n) * d) * cd + exp(1i * q(n) * d) * cd.';
    end
    Ghat(:, :, n) = G;
end
end
