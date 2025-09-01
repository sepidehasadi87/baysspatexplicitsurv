# BaySSpaTexplicitSurv

**BaySSpatexplicitSurv** is an R package for simulating and modeling spatial survival data. The package implements **Metropolis-Hastings MCMC** methods for parameter estimation in **STSU Weibull** models and provides tools for convergence diagnostics and model evaluation metrics such as DIC and BIC.

## Features

- Simulate spatial data with covariates and time-to-event outcomes.  
- Define priors and likelihoods for Weibull models with spatial structure.  
- Run **Metropolis-Hastings** algorithm for parameter estimation.  
- Compute model evaluation metrics including DIC, BIC, and MSE.  
- MCMC convergence diagnostics: Gelman-Rubin, trace plots, cumuplot, and correlation plots.

## Dependencies

The package requires the following R packages:

```r
install.packages(c(
  "MHadaptive", "MLmetrics", "metropolis", "coda", "BayesianTools",
  "geoR", "purrr", "DWreg", "DiscreteWeibull", "invgamma",
  "mvtnorm", "FAmle", "splines", "R2WinBUGS", "mcmcplots",
  "MCMCpack", "psych"
))
