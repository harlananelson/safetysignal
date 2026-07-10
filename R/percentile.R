# Posterior percentiles ----

#' Compute percentiles of 2-component Gamma mixture posterior
#'
#' For each drug-event pair, finds the value x such that the posterior
#' mixture CDF equals the given percentile. Uses [stats::uniroot()] for
#' root-finding.
#'
#' @param posterior_data A tibble from [compute_posterior()].
#' @param percentile Numeric in (0, 1). Default `0.05` for EB05.
#' @param col_name Name of the output column. Default derived from
#'   percentile (e.g., `"eb05"` for 0.05).
#'
#' @return The input tibble with an additional column for the percentile value.
#'
#' @export
posterior_percentile <- function(posterior_data, percentile = 0.05, col_name = NULL) {
  if (percentile <= 0 || percentile >= 1) {
    cli::cli_abort("{.arg percentile} must be between 0 and 1 (exclusive).")
  }

  if (is.null(col_name)) {
    col_name <- paste0("eb", sprintf("%02d", round(percentile * 100)))
  }

  idx <- seq_len(nrow(posterior_data))
  quantile_one <- function(i) {
    .mixture_quantile(
      p = percentile,
      q = posterior_data$q_post[i],
      alpha1 = posterior_data$alpha1_post[i],
      beta1 = posterior_data$beta1_post[i],
      alpha2 = posterior_data$alpha2_post[i],
      beta2 = posterior_data$beta2_post[i]
    )
  }

  # The per-pair root-find is the compute bottleneck and is embarrassingly
  # parallel: each pair is independent and .mixture_quantile is deterministic,
  # so forking across cores gives IDENTICAL results, just faster. Serial for
  # small inputs (fork overhead) or when cores == 1 (opt-out / non-unix).
  cores <- .ss_cores()
  vals <- if (cores > 1L && length(idx) >= .SS_PARALLEL_MIN) {
    res <- parallel::mclapply(idx, quantile_one, mc.cores = cores)
    vapply(
      res,
      function(x) if (length(x) == 1L && is.numeric(x)) as.double(x) else NA_real_,
      double(1)
    )
  } else {
    vapply(idx, quantile_one, double(1))
  }

  posterior_data[[col_name]] <- vals
  posterior_data
}

# Minimum number of pairs before parallelizing (fork overhead not worth it below).
.SS_PARALLEL_MIN <- 2000L

#' Resolve the worker-core count for the per-pair posterior quantile solve.
#'
#' Controlled by `options(safetysignal.cores = N)`; defaults to one fewer than
#' the detected cores. Forced to 1 on non-unix platforms (mclapply forking is
#' unix-only). Set to 1 to disable parallelism entirely.
#' @noRd
.ss_cores <- function() {
  n <- getOption("safetysignal.cores", NULL)
  if (is.null(n)) {
    n <- tryCatch(parallel::detectCores(), error = function(e) 1L)
    n <- max(1L, n - 1L)
  }
  if (.Platform$OS.type != "unix") n <- 1L
  max(1L, as.integer(n))
}


# Quantile of 2-component Gamma mixture via root-finding ----
.mixture_quantile <- function(p, q, alpha1, beta1, alpha2, beta2) {
  # Mixture CDF
  mix_cdf <- function(x) {
    q * pgamma(x, shape = alpha1, rate = beta1) +
      (1 - q) * pgamma(x, shape = alpha2, rate = beta2)
  }

  # Find a reasonable upper bound
  upper <- max(
    qgamma(0.9999, shape = alpha1, rate = beta1),
    qgamma(0.9999, shape = alpha2, rate = beta2)
  )
  # Ensure upper is large enough
  while (mix_cdf(upper) < p + 0.01) {
    upper <- upper * 2
  }

  result <- tryCatch(
    uniroot(
      \(x) mix_cdf(x) - p,
      lower = 1e-12,
      upper = upper,
      tol = 1e-8
    ),
    error = \(e) NULL
  )

  if (is.null(result)) {
    NA_real_
  } else {
    result$root
  }
}
