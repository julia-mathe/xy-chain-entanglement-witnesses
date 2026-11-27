% CHECK_QUENCH_TIME_EVOLUTION  omega_CMC(t) after the quench (Sec. 5, Fig. 7(a)).
% Needs CVX. Several minutes. Expected: omega_CMC crosses 1 near t = 0.19 and is about 0.83 at late times.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
require_cvx();

gamma = 1.0;  h_i = 0.3;  h_f = 2.0;  T_i = 1.0;
M0 = majorana_state_symbol_tdl('thermal', gamma, h_i, T_i);
times = [0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30];

for t = times
    tic;
    R = floor(min(30 + 5*t, 200));
    c = spin_correlations_symbol_tdl(majorana_state_symbol_tdl('evolved', M0, gamma, h_f, t), R);
    tail = max(max(max(abs(c(end-4:end, :, :)))));
    omega_sdp = cmc_block_sdp_tdl(c, 64);
    omega_formula = omega_structure_factors_tdl(c);
    fprintf('t = %5.2f: omega_SDP = %.4f  (structure factors %.4f), max|xy| = %.3f, tail = %.1e, R = %d  [%.0f s]\n', ...
        t, omega_sdp, omega_formula, max(abs(c(:, 1, 2))), tail, R, toc);
end
