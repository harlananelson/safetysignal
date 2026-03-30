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

  vals <- vapply(
    seq_len(nrow(posterior_data)),
    \(i) {
      .mixture_quantile(
        p = percentile,
        q = posterior_data$q_post[i],
        alpha1 = posterior_data$alpha1_post[i],
        beta1 = posterior_data$beta1_post[i],
        alpha2 = posterior_data$alpha2_post[i],
        beta2 = posterior_data$beta2_post[i]
      )
    },
    double(1)
  )

  posterior_data[[col_name]] <- vals
  posterior_data
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
