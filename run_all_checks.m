% RUN_ALL_CHECKS  Runs the check scripts one after the other; include_slow = true adds the slow ones.

include_slow = false;
here = fileparts(mfilename('fullpath'));
fast = {'check_thermodynamic_limit', 'check_theorem1_N10', 'check_theorem1_random_parameters', ...
    'check_free_fermions_open_chain', 'check_fully_ppt_example_N4', 'check_quench_temperature_window', ...
    'check_quench_reduced_states', 'check_quench_steady_state'};
slow = {'check_moment_hierarchy_N4', 'check_geometric_measure_certified', ...
    'check_geometric_measure_product_states', 'check_quench_time_evolution', 'check_quench_open_chain_N8', ...
    'check_quench_windows_and_open_chains', 'check_fully_ppt_open_chain_N4', 'check_quench_scan_Ti1'};
scripts = fast;
if include_slow, scripts = [fast, slow]; end
for n = 1:numel(scripts)
    run_one(fullfile(here, [scripts{n} '.m']));
end

function run_one(file)
[~, name] = fileparts(file);
fprintf('\n===== %s =====\n', name);
run(file);
end
