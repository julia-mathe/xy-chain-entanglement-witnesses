% CHECK_QUENCH_TEMPERATURE_WINDOW  Initial temperatures for which the quench generates entanglement
% (Sec. 5, Fig. 7(b)).
% No CVX needed. About a minute. Expected: steady state detected for T_i < 1.374, initial state for T_i < 0.718.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

gamma = 1.0;  h_i = 0.3;  h_f = 2.0;
omega_initial = @(T) omega_structure_factors_tdl(spin_correlations_symbol_tdl( ...
    majorana_state_symbol_tdl('thermal', gamma, h_i, T), 60));
omega_steady  = @(T) omega_structure_factors_tdl(spin_correlations_symbol_tdl( ...
    majorana_state_symbol_tdl('steady', majorana_state_symbol_tdl('thermal', gamma, h_i, T), gamma, h_f), 50));

options = optimset('TolX', 1e-5);
T_steady  = fzero(@(T) omega_steady(T) - 1, [1.30, 1.45], options);
T_initial = fzero(@(T) omega_initial(T) - 1, [0.68, 0.75], options);
fprintf('steady state (GGE) detected for T_i < %.3f;  initial state detected for T_i < %.3f\n', T_steady, T_initial);
