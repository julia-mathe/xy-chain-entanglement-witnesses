# Structure-factor witnesses of thermal entanglement in the XY chain — code

MATLAB code for

J. Mathé, A. Usui, O. Gühne and G. Vitagliano, *Structure-factor witnesses of thermal entanglement: exact optimal covariance criteria for the XY chain*, arXiv:2611.33284.

Section, figure, table and equation numbers below refer to this paper.

## Requirements

- MATLAB R2016b or newer; the figure scripts need R2020b or newer (`tiledlayout`, `exportgraphics`).
  The checks marked "no" in the CVX column also run in GNU Octave 8.
- [CVX](https://github.com/cvxr/CVX) for all semidefinite programs. The free solvers SDPT3 and SeDuMi
  that come with CVX suffice for the checks; MOSEK is recommended for the scans behind Figs. 2, 8 and 9.
- Fig. 13 only: the function `fpptwit` from the PPT-mixer code of Jungnitsch, Moroder and Gühne,
  Phys. Rev. Lett. 106, 190502 (2011), and YALMIP. It is not included here; without it the PPT mixer is skipped.

## Layout

| Folder | Content |
|---|---|
| `functions/` | building blocks (exact diagonalization, CMC, free fermions, quench, moment hierarchy, geometric measure) |
| `checks/` | scripts that reproduce the numbers quoted in the paper; each header lists the expected output |
| `figures/` | scripts that compute the data of the figures and draw them |
| `data/` | data files read by the figure scripts |
| `output/` | created on the first run; figures and large data files are written here |

Every script finds the other folders by itself; run it from any folder, e.g. `run checks/check_theorem1_N10.m`.
`checks/run_all_checks` runs the fast checks (set `include_slow = true` for all).

## Checks

| Script | Verifies | Paper | CVX |
|---|---|---|---|
| `check_theorem1_N10` | SDP = Theorem 1 at the XX, Ising and (0.5, 0.5) points, N = 10, PBC and OBC | Sec. 4.2, Tables 1–2, Fig. 4 | yes |
| `check_theorem1_random_parameters` | SDP = Theorem 1 at random (γ, h, T), N = 6 | Sec. 4.2 | yes |
| `check_thermodynamic_limit` | ω_CMC of the infinite chain | Table 1, Sec. 4.3, App. A | no |
| `check_free_fermions_open_chain` | free fermions vs. exact diagonalization; window of 10 in 200 spins | App. A, Sec. 4.1 | part 2 |
| `check_fully_ppt_example_N4` | fully PPT state detected by the CMC | Sec. 3, App. D | yes |
| `check_fully_ppt_open_chain_N4` | CMC vs. PPT (all cuts) in the (γ, h) plane, N = 4 | Sec. 3, App. D | yes |
| `check_moment_hierarchy_N4` | first level of the moment hierarchy = CMC for one- and two-spin data | App. C | yes |
| `check_geometric_measure_product_states` | Néel-type states maximize the inner problem of Eq. (36) | App. H | no |
| `check_geometric_measure_certified` | bounds on the geometric measure | Table 3, App. H | yes |
| `check_quench_steady_state` | initial state, steady state (GGE), equal-energy Gibbs state | Sec. 5, App. G | block SDP |
| `check_quench_temperature_window` | initial temperatures with generated entanglement | Sec. 5, Fig. 7(b) | no |
| `check_quench_time_evolution` | ω_CMC(t) | Sec. 5, Fig. 7(a) | yes |
| `check_quench_reduced_states` | partial transposes of two- and three-spin blocks | Fig. 14(b) | no |
| `check_quench_open_chain_N8` | open chain of 8 spins: all cuts and the CMC | App. G | part 2 |
| `check_quench_windows_and_open_chains` | windows of L spins and open chains | Fig. 14(c) | yes |
| `check_quench_scan_Ti1` | which quenches generate entanglement | Sec. 5, Fig. 14(a) | no |

## Figures

| Figure | Scripts | Data in `data/` |
|---|---|---|
| 1 | schematic; the witness types follow from `check_theorem1_N10` | – |
| 2, 8, 13 | `fig02_08_13_data` (CVX, MOSEK), then `fig02_08_13_criteria_vs_temperature` | `N4_criteria_vs_T_<BC>.mat` |
| 3 | `fig03_witness_type_tdl` | – |
| 4, 11, 12 | `fig04_11_12_data` (CVX), then `fig04_11_12_matrices` | `matrix_figures_N10_T0.1.mat` |
| 5 | `fig05_structure_factors` | – |
| 6 | `fig06_data_*`, then `fig06_entanglement_temperature` | `entanglement_temperature_*.mat`, `energy_witness_threshold.mat`, `omega_vs_temperature.mat` |
| 7, 14 | `fig07_14_quench` (CVX) | – |
| 9 | `fig09_phase_diagram_N4` (CVX, MOSEK; once per boundary condition) | – |
| 10 | `fig10_cmc_vs_ppt_N12` (CVX) | – |

## Conventions

J = 1 and k_B = 1. Spins are numbered 1, …, N, spin 1 being the leftmost tensor factor. In covariance
matrices σ_j^α has index 3(j−1)+α with α = 1, 2, 3 for x, y, z. Ring momenta are q = 2πm/N.
