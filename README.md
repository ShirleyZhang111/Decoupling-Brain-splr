# Decoupling whole-brain dynamics (Sparse-plus-Low-Rank Effective Connectivity)

This repository contains MATLAB code accompanying a **complex-valued generative framework** that decomposes whole-brain *effective coupling* into:

- a **low-rank (global) component** capturing coherent, brain-wide fluctuations, and  
- a **sparse component** capturing functionally specific interactions.

The implementation supports two coupling mechanisms:

- **`bipartite`**: bipartite conjugate coupling (Hamiltonian-inspired model)  
- **`amplitude`**: amplitude-based coupling (real-valued baseline, traditional Hopf model)

## Repository structure

- `demo_sparse_plus_low_rank_model.m` — minimal demo for fitting the model on example rs-fMRI data.
- `generate_coupling_components.m` — main entry point for estimating coefficients and coupling components.
- `construct_data.m` — constructs model-predicted derivatives from polynomial terms and coupling.
- `calc_GSR_R2.m` — computes variance explained by global signal regression (R² per parcel/region).
- `analysis/` — scripts used for group-level and task analyses:
  - `analysis_group_level_rest.m`
  - `analysis_group_level_tasks.m`
  - `analysis_language_coupling_circuits.m`
- `functions/` — optimization, metrics, and utilities:
  - ADMM solvers: `hopf_admm_*`
  - graph/cycle utilities: `top_cycles_exactL_sumabs.m`, etc.
  - network metrics: `calc_network_metrics.m`, `calc_node_wise_metrics.m`
- `data/` — small example `.mat` files used by the demo and analysis scripts.

## Requirements

- MATLAB (tested with recent versions).
- Sufficient RAM for loading `data/*.mat` and fitting ADMM models on your parcellation size.

## Quick start

1. **Clone / download** this repository.
2. Start MATLAB and `cd` into the repository root.
3. Run the demo:

```matlab
demo_sparse_plus_low_rank_model
```

By default, the demo loads:

- `./data/HCPex_rest_data.mat` (expects a variable `TC`, shaped `[N x T]`)

and writes results to:

- `./results/` (created if not present)

### Key configuration options (from the demo)

```matlab
cfg.rank = 1;              % rank of low-rank component (0,1,2,...)
cfg.maxit = 2000;          % ADMM iterations
cfg.tol  = 1e-12;          % stopping tolerance
cfg.mechanism = 'bipartite'; % 'bipartite' or 'amplitude'

cfg.T  = 101:600;          % time indices used for fitting
cfg.mu = 0.1*length(cfg.T);% sparsity penalty (larger => sparser)
cfg.xi = 100;              % low-rank regularization weight (see code)

cfg.data_path = './data/';
cfg.file_name = 'HCPex_rest_data.mat';
cfg.save_path = './results/';
```

The main fitting call is:

```matlab
[coeff, sparse, u, v, v_hat, res, mse_real, cos_theta] = generate_coupling_components(cfg);
```

### Outputs (typical)

Depending on the chosen `rank` and `mechanism`, the main function returns:

- `coeff` — coefficients of the polynomial basis in the local dynamics
- `sparse` — estimated sparse coupling matrix
- `u, v` — low-rank factors (so the low-rank coupling is approximately `u*v'`)
- `res` — residuals
- `mse_real` — mean squared error of the data fit
- `cos_theta` — subspace similarity metric (for stability checks)

## Data format

The example datasets in `data/` are MATLAB `.mat` files that include a variable:

- `TC`: time series matrix of shape **`[N x T]`**, where `N` is number of parcels/regions and `T` is number of time points.

The demo trims the first/last few time points internally (see `calc_GSR_R2.m` and related scripts), so ensure your data have adequate length.

## Reproducing analyses

The scripts in `analysis/` are intended to reproduce the main figures/analyses reported in the paper. Typical usage is:

1. Fit individual-level models and save outputs
2. Run group-level aggregation scripts in `analysis/`

Because paths and filenames are dataset-specific, you may need to edit:

- `data_path`, `file_name`, and output directories in each analysis script.

## Notes on interpretation

- The inferred matrices represent **BOLD-level directed dynamical dependencies**, not synaptic or causal ground truth.
- The **low-rank component** is intended to capture a dominant global coordination mode.
- The **sparse component** captures residual, network-specific interactions after accounting for the global background.

## Citation

If you use this code, please cite the accompanying manuscript (add your BibTeX entry here).

## License

See [LICENSE](LICENSE).
