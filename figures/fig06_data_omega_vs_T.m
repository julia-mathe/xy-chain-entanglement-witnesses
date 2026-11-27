% FIG06_DATA_OMEGA_VS_T  omega_CMC versus T at the XX and Ising points (Fig. 6(b)).
% No CVX needed.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

temperatures = [linspace(0.05, 0.5, 10), linspace(0.6, 2.6, 21)];
points = {0, 0, 'XX'
          1, 1, 'Ising'};

omega_infinite = zeros(2, numel(temperatures));
omega_N10 = zeros(2, numel(temperatures));
for n = 1:2
    gamma = points{n, 1};
    h = points{n, 2};
    ring = ring_spectrum(10, gamma, h);
    for it = 1:numel(temperatures)
        omega_infinite(n, it) = cmc_tdl(gamma, h, temperatures(it));
        omega_N10(n, it) = ring_cmc_at_temperature(ring, temperatures(it));
    end
    fprintf('%s point, infinite chain, first 12 temperatures: %s\n', points{n, 3}, ...
            sprintf('%.3f ', omega_infinite(n, 1:12)));
end

save(fullfile(repo, 'data', 'omega_vs_temperature.mat'), 'temperatures', 'omega_infinite', 'omega_N10', '-v7');
fprintf('Saved data/omega_vs_temperature.mat\n');
