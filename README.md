# Habitat suitability, step selection, and connectivity of Asiatic black bears on Shikoku

R scripts supporting the landscape-scale resource selection function (RSF), integrated step-selection function (iSSF), model validation, and connectivity simulations in a *Biological Conservation* manuscript.

## Contents

- **RSF (`R/01–04`):** Landscape preparation, model fitting, suitability prediction, and Figures 4–5.
- **iSSF (`R/10–14`):** Model fitting, leave-one-individual-out validation, Figure 6, and raster-resolution sensitivity analysis.
- **Connectivity (`R/20–21`):** iSSF-based movement simulations and connectivity summaries.
- **`R/00_config.R`:** Input and output paths shared across scripts.

## Usage

Run the scripts from the repository root using R. Required input files should be placed under `data_private/` using the paths specified in `R/00_config.R`. The iSSF scripts use the preprocessed 100-m step dataset.

```bash
# Landscape RSF
code/RSF/01_prepare_landscape.R
code/RSF/02_fit_predict_rsf.R
code/RSF/03_figure4_coefficients.R
code/RSF/04_figure5_covariates.R

# iSSF
code/iSSF/11_fit_issf.R
code/iSSF/12_validate_issf.R
code/iSSF/13_figure6_logrss.R
code/iSSF/14_resolution_sensitivity.R  # optional

# Connectivity simulations
code/Simulation/21_run_connectivity.R --smoke-only
code/Simulation/21_run_connectivity.R
```

Scripts `10_issf_common.R` and `20_simulation_functions.R` contain functions called by other scripts. Generated files are saved under `results/`.

## Data availability

Bear GPS locations and associated spatial datasets are not publicly available because of conservation sensitivity and data-sharing restrictions. This repository contains analysis code only; restricted data and generated spatial outputs are excluded.
