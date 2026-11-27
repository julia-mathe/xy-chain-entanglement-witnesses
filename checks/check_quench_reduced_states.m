% CHECK_QUENCH_REDUCED_STATES  Partial transposes of two- and three-spin reduced states after the quench
% (Sec. 5, Fig. 14(b)).
% No CVX needed. Expected: all pairs PPT; contiguous three-spin blocks NPT from t = 0.2 on and in the GGE.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

gamma = 1.0;  h_i = 0.3;  h_f = 2.0;  T_i = 1.0;
M0 = majorana_state_symbol_tdl('thermal', gamma, h_i, T_i);
times = [0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1, 2, 3, 5, 10, 20];

for n = 1:numel(times) + 1
    if n <= numel(times)
        Mh = majorana_state_symbol_tdl('evolved', M0, gamma, h_f, times(n));
        label = sprintf('t = %5.2f', times(n));
    else
        Mh = majorana_state_symbol_tdl('steady', M0, gamma, h_f);
        label = 'GGE     ';
    end
    Mw = majorana_window_tdl(Mh, 4);
    pair_min = zeros(1, 3);
    for d = 1:3
        rho2 = reduced_density_matrix_gaussian(Mw, [1, 1 + d]);
        pair_min(d) = min(real(eig(partial_transpose(rho2, 1, 2))));
    end
    rho3 = reduced_density_matrix_gaussian(Mw, [1, 2, 3]);
    block_min = Inf;
    for A = 1:3
        block_min = min(block_min, min(real(eig(partial_transpose(rho3, A, 3)))));
    end
    fprintf('%s: smallest PT eigenvalue of pairs d = 1, 2, 3: [%7.4f %7.4f %7.4f]   three-spin block: %+.4f\n', ...
        label, pair_min, block_min);
end
