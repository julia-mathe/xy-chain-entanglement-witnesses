% FIG06_DATA_ENERGY_WITNESS  Detection temperature of the energy witness (Fig. 6(a), Sec. 4.3).
% No CVX needed.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

fprintf('XX point T_E = %.5f,   Ising point T_E = %.5f\n', ...
        energy_witness_temperature(0, 0), energy_witness_temperature(1, 1));

fields = linspace(0, 2.5, 51);
gammas = [0, 1];
T_E = zeros(numel(gammas), numel(fields));
for ig = 1:numel(gammas)
    for ih = 1:numel(fields)
        T_E(ig, ih) = energy_witness_temperature(gammas(ig), fields(ih));
    end
    fprintf('gamma = %d, h = 0, 0.25, ..., 2.5:  %s\n', gammas(ig), sprintf('%.3f ', T_E(ig, 1:5:end)));
end

save(fullfile(repo, 'data', 'energy_witness_threshold.mat'), 'gammas', 'fields', 'T_E', '-v7');
fprintf('Saved data/energy_witness_threshold.mat\n');
