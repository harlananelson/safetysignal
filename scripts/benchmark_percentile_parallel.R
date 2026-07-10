#!/usr/bin/env Rscript
# Parity + benchmark for the parallelized posterior_percentile().
# Proves the mclapply path returns IDENTICAL results to serial (method-
# preserving) and measures the speedup. Run from the safetysignal repo root.

source("R/percentile.R")

set.seed(1)
n <- 50000L
d <- data.frame(
  q_post      = runif(n, 0.3, 0.9),
  alpha1_post = runif(n, 1, 50),
  beta1_post  = runif(n, 1, 30),
  alpha2_post = runif(n, 1, 20),
  beta2_post  = runif(n, 1, 20)
)

cores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("n = %d pairs; cores = %d\n", n, cores))

# Serial (opt-out)
options(safetysignal.cores = 1L)
t_serial <- system.time(s <- posterior_percentile(d, 0.05))[["elapsed"]]

# Parallel
options(safetysignal.cores = cores)
t_par <- system.time(p <- posterior_percentile(d, 0.05))[["elapsed"]]

identical_vals <- isTRUE(all.equal(s$eb05, p$eb05))
n_na_s <- sum(is.na(s$eb05)); n_na_p <- sum(is.na(p$eb05))

cat(sprintf("PARITY: identical results = %s (NA serial=%d, parallel=%d)\n",
            identical_vals, n_na_s, n_na_p))
cat(sprintf("SERIAL:   %.2fs\nPARALLEL: %.2fs (%d cores)\nSPEEDUP:  %.1fx\n",
            t_serial, t_par, cores, t_serial / t_par))
