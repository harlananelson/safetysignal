# safetysignal — Bayesian GPS Signal Detection

## Overview
R package implementing the Multi-item Gamma-Poisson Shrinker (MGPS) for pharmacovigilance signal detection. Uses a 2-component Gamma mixture prior fitted by empirical Bayes (marginal likelihood of the counts, DuMouchel 1999 / openEBGM), with full posterior percentiles (not EBGM approximation).

## Architecture

```
R/
├── gamma-mle.R          # Newton-Raphson for Gamma shape MLE
├── observed-expected.R  # O/E table from contingency data
├── prior.R              # Empirical Bayes fit of the 2-component Gamma mixture (+ calibration check)
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

- **Prior fit**: `fit_prior()` maximises the per-component zero-truncated negative-binomial marginal likelihood of (observed, expected), matching `openEBGM::autoHyper`. Never fit the prior to raw O/E ratios: in sparse data single reports with tiny E have ratios in the thousands and inflate EB05 for low-count pairs. The old ratio-EM fit survives only as `method = "rr_em"`. Use `check_prior_calibration()` to confirm the prior reproduces the count distribution.

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
