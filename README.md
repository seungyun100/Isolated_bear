# Habitat suitability, step selection, and connectivity of Asiatic black bears on Shikoku

R code supporting the landscape-scale resource selection function (RSF), integrated step-selection function (iSSF), and connectivity simulations described in a *Biological Conservation* manuscript.

## Repository structure

```text
code/
  01_config.R
  RSF/
    01_prepare_landscape.R
    02_fit_predict_rsf.R
    03_figure4_coefficients.R
    04_figure5_covariates.R
  iSSF/
    01_issf_common.R
    02_fit_issf.R
    03_validate_issf.R
    04_figure6_logrss.R
    05_resolution_sensitivity.R
  Simulation/
    01_simulation_functions.R
    02_run_connectivity.R
```

## Usage

Run the scripts from the repository root with R. Configure input paths in `code/01_config.R` and place authorized input data in `data_private/`. The iSSF analysis uses preprocessed 100-m step data.

```bash
# Landscape RSF
Rscript code/RSF/01_prepare_landscape.R
Rscript code/RSF/02_fit_predict_rsf.R
Rscript code/RSF/03_figure4_coefficients.R
Rscript code/RSF/04_figure5_covariates.R

# iSSF
Rscript code/iSSF/02_fit_issf.R
Rscript code/iSSF/03_validate_issf.R
Rscript code/iSSF/04_figure6_logrss.R
Rscript code/iSSF/05_resolution_sensitivity.R  # optional

# Connectivity simulations
Rscript code/Simulation/02_run_connectivity.R --smoke-only
Rscript code/Simulation/02_run_connectivity.R
```

## Data availability

Bear GPS locations and associated spatial data are not publicly available because of conservation sensitivity and data-sharing restrictions. This repository provides analysis code only.
