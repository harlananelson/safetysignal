# Posterior computation ----

#' Compute posterior Gamma mixture for each drug-event pair
#'
#' Given observed counts and expected counts from [compute_observed_expected()],
#' and a fitted prior from [fit_prior()], computes the 2-component Gamma
#' mixture posterior for each pair.
#'
#' The posterior for component k has:
#' - `alpha_k_post = alpha_k + observed`
#' - `beta_k_post = beta_k + expected`
#' - Weight updated via Negative Binomial marginal likelihood
#'
#' @param oe_data A tibble from [compute_observed_expected()].
#' @param prior A `gps_prior` object from [fit_prior()].
#'
#' @return The input tibble with additional columns:
#'   `alpha1_post`, `beta1_post`, `alpha2_post`, `beta2_post`,
#'   `q_post` (posterior weight on component 1).
#'
#' @export
compute_posterior <- function(oe_data, prior) {
  if (!inherits(prior, "gps_prior")) {
    cli::cli_abort("{.arg prior} must be a {.cls gps_prior} object from {.fn fit_prior}.")
  }

  oe_data |>
    dplyr::mutate(
      # Posterior shape and rate
      alpha1_post = prior$alpha1 + .data$observed,
      beta1_post = prior$beta1 + .data$expected,
      alpha2_post = prior$alpha2 + .data$observed,
      beta2_post = prior$beta2 + .data$expected,

      # Posterior mixing weights via NB marginal likelihood
      # P(n | component k) = NB(n, size = alpha_k, prob = beta_k / (beta_k + E))
      log_nb1 = dnbinom(
        .data$observed,
        size = prior$alpha1,
        prob = prior$beta1 / (prior$beta1 + .data$expected),
        log = TRUE
      ) + log(prior$pi),
      log_nb2 = dnbinom(
        .data$observed,
        size = prior$alpha2,
        prob = prior$beta2 / (prior$beta2 + .data$expected),
        log = TRUE
      ) + log(1 - prior$pi),

      # Normalize (log-sum-exp)
      log_nb_max = pmax(.data$log_nb1, .data$log_nb2),
      q_post = exp(.data$log_nb1 - .data$log_nb_max) /
        (exp(.data$log_nb1 - .data$log_nb_max) + exp(.data$log_nb2 - .data$log_nb_max))
    ) |>
    dplyr::select(
      -"log_nb1", -"log_nb2", -"log_nb_max"
    )
}
