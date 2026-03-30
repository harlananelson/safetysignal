# Full pipeline convenience function ----

#' Run full GPS signal detection pipeline
#'
#' Convenience wrapper that runs the complete Gamma-Poisson Shrinker pipeline:
#' observed/expected → fit prior → posterior → percentiles → signal detection.
#'
#' @param data A data frame of drug-event reports.
#' @param drug_col Column identifying the drug/vaccine (tidy-select).
#' @param event_col Column identifying the adverse event (tidy-select).
#' @param count_col Optional column with pre-aggregated counts.
#' @param threshold Signal detection threshold for EB05. Default `1`.
#' @param prior Optional pre-fitted `gps_prior` object. If `NULL`, fits
#'   the prior from the data.
#' @param prior_tol Convergence tolerance for prior fitting. Default `1e-6`.
#'
#' @return A tibble with all computed columns: `drug`, `event`, `observed`,
#'   `expected`, `rr`, posterior parameters, `eb05`, `eb95`, `is_signal`,
#'   `signal_strength`.
#'
#' @export
gps_detect <- function(data, drug_col, event_col, count_col = NULL,
                       threshold = 1, prior = NULL, prior_tol = 1e-6) {
  drug_col <- rlang::enquo(drug_col)
  event_col <- rlang::enquo(event_col)
  count_col <- rlang::enquo(count_col)

  cli::cli_h3("GPS Signal Detection Pipeline")

  # Step 1: Observed/Expected
  cli::cli_inform("Computing observed/expected table...")
  if (rlang::quo_is_null(count_col)) {
    oe <- compute_observed_expected(data, !!drug_col, !!event_col)
  } else {
    oe <- compute_observed_expected(data, !!drug_col, !!event_col, !!count_col)
  }
  cli::cli_inform("{nrow(oe)} drug-event pairs.")

  # Step 2: Fit prior
  if (is.null(prior)) {
    cli::cli_inform("Fitting 2-component Gamma mixture prior via EM...")
    prior <- fit_prior(oe, tol = prior_tol)
    print(prior)
  } else {
    cli::cli_inform("Using pre-fitted prior.")
  }

  # Step 3-5: Posterior, percentiles, signals
  cli::cli_inform("Computing posteriors and percentiles...")
  result <- oe |>
    compute_posterior(prior) |>
    posterior_percentile(percentile = 0.05) |>
    posterior_percentile(percentile = 0.95) |>
    detect_signals(threshold = threshold)

  n_signals <- sum(result$is_signal, na.rm = TRUE)
  cli::cli_inform("{n_signals} signals detected (EB05 >= {threshold}).")

  result
}
