% CHECK_QUENCH_STEADY_STATE  Quench h_i = 0.3 -> h_f = 2 at T_i = 1, infinite Ising chain (Sec. 5, App. G).
% Block SDP needs CVX. A few minutes. Expected: initial 1.0432 (s = 0.155 0.978 0.953);
% GGE 0.9019 (s = 0.428 0.686 0.691; block SDP 0.9019); T_eff = 4.054, Gibbs state 1.2126.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
quick = false;
has_cvx = exist('cvx_begin', 'file') > 0;

%% Part 0: thermal states (validation against Table 1)
if quick, R0 = 100; else, R0 = 160; end
cases = [1.0, 1.0, 0.1, 0.559
         0.0, 0.0, 0.1, 0.189];
for n = 1:size(cases, 1)
    Mh = majorana_state_symbol_tdl('thermal', cases(n,1), cases(n,2), cases(n,3));
    c = spin_correlations_symbol_tdl(Mh, R0);
    [s, qmin, omega] = structure_factor_minima_tdl(c);
    fprintf('thermal (gamma,h,T) = (%g,%g,%g): omega = %.4f (paper %.3f), s =[%.4f %.4f %.4f], q/pi = [%.3f %.3f %.3f], max|xy| = %.1e\n', ...
        cases(n,1:3), omega, cases(n,4), s, qmin/pi, max(abs(c(:,1,2))));
end

%% Part 1: the quench of Sec. 5
gamma = 1.0;  h_i = 0.3;  h_f = 2.0;  T_i = 1.0;  R = 60;
M0 = majorana_state_symbol_tdl('thermal', gamma, h_i, T_i);
c0 = spin_correlations_symbol_tdl(M0, R);
[s0, q0, omega0] = structure_factor_minima_tdl(c0);
fprintf('\ninitial: omega = %.4f, s = [%.3f %.3f %.3f], q/pi = [%.3f %.3f %.3f]\n', omega0, s0, q0/pi);

Mg = majorana_state_symbol_tdl('steady', M0, gamma, h_f);
cg = spin_correlations_symbol_tdl(Mg, R);
[sg, qg, omegag] = structure_factor_minima_tdl(cg);
fprintf('GGE: omega (formula) = %.4f, s = [%.3f %.3f %.3f], q/pi = [%.3f %.3f %.3f], max|xy corr| = %.1e\n', ...
    omegag, sg, qg/pi, max(abs(cg(:,1,2))));
if has_cvx
    fprintf('GGE: block SDP (nq = 64) = %.4f\n', cmc_block_sdp_tdl(cg, 64));
else
    fprintf('GGE: block SDP skipped (CVX not found)\n');
end

e0 = energy_density_symbol_tdl(M0, gamma, h_f);
fprintf('energy density after the quench (H_f): %.6f, in the GGE: %.6f\n', e0, energy_density_symbol_tdl(Mg, gamma, h_f));

energy_difference = @(T) energy_density_symbol_tdl(majorana_state_symbol_tdl('thermal', gamma, h_f, T), gamma, h_f) - e0;
T_eff = fzero(energy_difference, [0.05, 50]);
cT = spin_correlations_symbol_tdl(majorana_state_symbol_tdl('thermal', gamma, h_f, T_eff), R);
[sT, ~, omegaT] = structure_factor_minima_tdl(cT);
fprintf('T_eff = %.3f, omega(Gibbs state of H_f at T_eff) = %.4f, s = [%.3f %.3f %.3f]\n', T_eff, omegaT, sT);
