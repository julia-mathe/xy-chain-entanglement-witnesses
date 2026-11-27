% CHECK_QUENCH_SCAN_TI1  Which quenches (h_i, h_f) at T_i = 1 generate entanglement (Sec. 5, Fig. 14(a)).
% No CVX needed. 10-30 min. Writes output/quench_scan_Ti1.mat.
% Expected: for h_i <= 0.4 only for h_f >= 1.4; at h_i = 0.5 also for h_f = 0.6 ... 1.0.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
if ~exist(fullfile(repo, 'output'), 'dir'), mkdir(fullfile(repo, 'output')); end

gamma = 1.0;  T_i = 1.0;  R = 40;
h_i_list = 0:0.1:3;
h_f_list = 0:0.1:3;
omega_initial = zeros(numel(h_i_list), 1);
omega_steady = zeros(numel(h_i_list), numel(h_f_list));

tic;
for a = 1:numel(h_i_list)
    M0 = majorana_state_symbol_tdl('thermal', gamma, h_i_list(a), T_i);
    omega_initial(a) = omega_structure_factors_tdl(spin_correlations_symbol_tdl(M0, R));
    for b = 1:numel(h_f_list)
        Mg = majorana_state_symbol_tdl('steady', M0, gamma, h_f_list(b));
        omega_steady(a, b) = omega_structure_factors_tdl(spin_correlations_symbol_tdl(Mg, R));
    end
    fprintf('h_i = %.1f done [%.0f s]\n', h_i_list(a), toc);
end
save(fullfile(repo, 'output', 'quench_scan_Ti1.mat'), 'h_i_list', 'h_f_list', 'omega_initial', ...
     'omega_steady', 'gamma', 'T_i', 'R', '-v7');

%% Summary
fprintf('\ninitial state detected for h_i in: %s\n', mat2str(h_i_list(omega_initial < 1)));
generated = (omega_steady < 1) & (omega_initial >= 1);
fprintf('generated (steady state detected, initial state not):\n');
for a = 1:numel(h_i_list)
    hf = h_f_list(generated(a, :));
    if ~isempty(hf)
        fprintf('  h_i = %.1f: h_f in [%.1f, %.1f] (%d points)\n', h_i_list(a), min(hf), max(hf), numel(hf));
    end
end
lost = (omega_steady >= 1) & (omega_initial < 1);
fprintf('lost (initial state detected, steady state not): %d of %d quenches\n', nnz(lost), numel(lost));
for a = 1:numel(h_i_list)
    hf = h_f_list(lost(a, :));
    if ~isempty(hf)
        fprintf('  h_i = %.1f: h_f = %s\n', h_i_list(a), mat2str(hf));
    end
end
