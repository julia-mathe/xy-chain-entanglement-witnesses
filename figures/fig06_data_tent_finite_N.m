% FIG06_DATA_TENT_FINITE_N  T_ent of rings with N = 4 and N = 10 along the XX and Ising lines (Fig. 6(a)).
% No CVX needed. A few minutes.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

fields = linspace(0, 2.5, 51);
gammas = [0, 1];

for N = [4 10]
    T_ent = zeros(numel(gammas), numel(fields));
    switches = zeros(numel(gammas), numel(fields));
    smallest_omega = zeros(numel(gammas), numel(fields));
    tic;
    for ig = 1:numel(gammas)
        for ih = 1:numel(fields)
            ring = ring_spectrum(N, gammas(ig), fields(ih));
            omega_of_T = @(T) ring_cmc_at_temperature(ring, T);
            [T_ent(ig, ih), temperatures, values] = entanglement_temperature_scan(omega_of_T, 0.02, 5, 40);
            detected = values < 1 - 1e-12;
            switches(ig, ih) = sum(diff(detected) ~= 0);
            smallest_omega(ig, ih) = min(values);
            if smallest_omega(ig, ih) > 1 - 1e-12
                T_ent(ig, ih) = 0;
            end
        end
        fprintf('N = %d, gamma = %d: T_ent(h = 0) = %.4f, T_ent(h = 1) = %.4f  (%.0f s)\n', ...
                N, gammas(ig), T_ent(ig, 1), T_ent(ig, fields == 1), toc);
    end
    file_name = fullfile(repo, 'data', sprintf('entanglement_temperature_N%d.mat', N));
    save(file_name, 'N', 'gammas', 'fields', 'T_ent', 'switches', 'smallest_omega', '-v7');
    fprintf('Saved %s\n', file_name);
end
