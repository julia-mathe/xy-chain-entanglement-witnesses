% FIG06_DATA_TENT_TDL  T_ent of the infinite chain along the XX and Ising lines (Fig. 6(a), Sec. 4.3).
% No CVX needed. A few minutes.

%% Setup
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));

fields = linspace(0, 2.5, 51);
gammas = [0, 1];
T_ent = zeros(numel(gammas), numel(fields));
momenta_at_T_ent = nan(numel(gammas), numel(fields), 3);

tic;
for ig = 1:numel(gammas)
    for ih = 1:numel(fields)
        gamma = gammas(ig);
        h = fields(ih);
        omega_of_T = @(T) cmc_tdl(gamma, h, T);
        T_ent(ig, ih) = entanglement_temperature_from_above(omega_of_T, 4.0, 0.25, 24);
        if T_ent(ig, ih) > 0
            [~, ~, qmin] = cmc_tdl(gamma, h, T_ent(ig, ih));
            momenta_at_T_ent(ig, ih, :) = qmin;
        end
        fprintf('gamma = %d, h = %.2f: T_ent = %.4f  (%.0f s)\n', gamma, h, T_ent(ig, ih), toc);
    end
end

save(fullfile(repo, 'data', 'entanglement_temperature_tdl.mat'), 'gammas', 'fields', 'T_ent', ...
     'momenta_at_T_ent', '-v7');
fprintf('Saved data/entanglement_temperature_tdl.mat\n');

small = fields >= 0.05 & fields <= 0.2;
fprintf('Ising line, 0.05 <= h <= 0.2:  T_ent * ln(1/h) = %s\n', ...
        sprintf('%.3f ', T_ent(2, small) .* log(1 ./ fields(small))));
