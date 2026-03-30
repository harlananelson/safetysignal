# safetysignal — Bayesian GPS Signal Detection

## Overview
R package implementing the Multi-item Gamma-Poisson Shrinker (MGPS) for pharmacovigilance signal detection. Uses a 2-component Gamma mixture prior fitted via EM, with full posterior percentiles (not EBGM approximation).

## Architecture

```
R/
├── gamma-mle.R          # Newton-Raphson for Gamma shape MLE
├── observed-expected.R  # O/E table from contingency data
├── prior.R              # EM for 2-component Gamma mixture
├── posterior.R           # Posterior weights via NB marginal likelihood
├── percentile.R          # Mixture quantiles via uniroot
├── signal.R              # Signal classification
└── gps-detect.R          # Full pipeline convenience function
```

## Key Design Decisions

- **Full posterior, not EBGM**: Signal detection uses percentiles of the 2-component Gamma mixture posterior directly. The EBGM geometric mean is an unnecessary approximation.
- **Gamma parameterization**: Uses `shape, rate` throughout (not `shape, scale`), matching DuMouchel's notation.
- **Component ordering**: Component 1 is always "background" (smaller mean), component 2 is "signal" (larger mean).

## Downstream Apps

- **faers-mobi** — Vaccine safety (VAERS data)
- **aers-mobi** — Drug/device safety (FAERS + MAUDE data)
- **globalpatientsafety.com** — Broader safety metrics

## Statistical Method

Prior: π·Gamma(α₁,β₁) + (1-π)·Gamma(α₂,β₂)
Posterior component k: Gamma(αₖ + O, βₖ + E)
Posterior weights: via NB(O | αₖ, βₖ/(βₖ+E)) marginal likelihood
Signal: EB05 (5th percentile of posterior) > threshold

## Conventions

- R >= 4.1, native pipe |>
- tidyverse style (snake_case, <- assignment)
- testthat 3e
- cli for user-facing messages
- roxygen2 for documentation
