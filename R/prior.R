# EM algorithm for 2-component Gamma mixture prior ----

#' Fit 2-component Gamma mixture prior via EM
#'
#' Fits `pi * Gamma(alpha1, beta1) + (1 - pi) * Gamma(alpha2, beta2)`
#' to the observed RR values using Expectation-Maximization.
#'
#' @param data A tibble from [compute_observed_expected()] (uses `rr` column),
#'   or a numeric vector of observed/expected ratios.
#' @param tol Convergence tolerance on log-likelihood. Default `1e-6`.
#' @param max_iter Maximum EM iterations. Default 200.
#' @param init Optional named list of initial parameters:
#'   `list(alpha1, beta1, alpha2, beta2, pi)`. If `NULL`, uses
#'   moment-based initialization.
#'
#' @return A list of class `"gps_prior"` with components:
#'   `alpha1`, `beta1`, `alpha2`, `beta2`, `pi`, `log_lik`,
#'   `converged`, `iterations`.
#'
#' @export
fit_prior <- function(data, tol = 1e-6, max_iter = 200, init = NULL) {
  # Extract RR values
  if (is.data.frame(data)) {
    if (!"rr" %in% names(data)) {
      cli::cli_abort("{.arg data} must have a column named {.field rr}.")
    }
    rr <- data$rr
  } else {
    rr <- as.numeric(data)
  }

  # Remove zeros and NAs
  rr <- rr[is.finite(rr) & rr > 0]
  n <- length(rr)

  if (n < 10) {
    cli::cli_abort("Need at least 10 positive RR values to fit prior, got {n}.")
  }

  # Initialize parameters
  if (is.null(init)) {
    init <- .init_prior(rr)
  }

  alpha1 <- init$alpha1
  beta1 <- init$beta1
  alpha2 <- init$alpha2
  beta2 <- init$beta2
  pi_mix <- init$pi

  log_lik_prev <- -Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    # E-step: posterior responsibilities
    log_d1 <- dgamma(rr, shape = alpha1, rate = beta1, log = TRUE) + log(pi_mix)
    log_d2 <- dgamma(rr, shape = alpha2, rate = beta2, log = TRUE) + log(1 - pi_mix)

    # Log-sum-exp for numerical stability
    log_dmax <- pmax(log_d1, log_d2)
    log_sum <- log_dmax + log(exp(log_d1 - log_dmax) + exp(log_d2 - log_dmax))

    r <- exp(log_d1 - log_sum)  # responsibility for component 1

    # Log-likelihood
    log_lik <- sum(log_sum)

    # Check convergence
    if (abs(log_lik - log_lik_prev) < tol) {
      converged <- TRUE
      break
    }
    log_lik_prev <- log_lik

    # M-step
    pi_mix <- mean(r)

    # Clamp pi away from 0/1 to prevent degenerate solutions
    pi_mix <- max(0.01, min(0.99, pi_mix))

    # Fit each component via weighted Gamma MLE
    fit1 <- gamma_shape_mle(rr, weights = r)
    fit2 <- gamma_shape_mle(rr, weights = 1 - r)

    alpha1 <- fit1$alpha
    beta1 <- fit1$beta
    alpha2 <- fit2$alpha
    beta2 <- fit2$beta
  }

  # Ensure component 1 is the "background" (higher rate = smaller mean)
  if (alpha1 / beta1 > alpha2 / beta2) {
    result <- list(
      alpha1 = alpha2, beta1 = beta2,
      alpha2 = alpha1, beta2 = beta1,
      pi = 1 - pi_mix
    )
  } else {
    result <- list(
      alpha1 = alpha1, beta1 = beta1,
      alpha2 = alpha2, beta2 = beta2,
      pi = pi_mix
    )
  }

  result$log_lik <- log_lik
  result$converged <- converged
  result$iterations <- iter
  result$n <- n
  class(result) <- "gps_prior"
  result
}


#' @export
print.gps_prior <- function(x, ...) {
  cli::cli_h3("GPS 2-Component Gamma Mixture Prior")
  cli::cli_bullets(c(
    "*" = "Component 1 (background): Gamma({round(x$alpha1, 4)}, {round(x$beta1, 4)}) — mean = {round(x$alpha1 / x$beta1, 4)}",
    "*" = "Component 2 (signal):     Gamma({round(x$alpha2, 4)}, {round(x$beta2, 4)}) — mean = {round(x$alpha2 / x$beta2, 4)}",
    "*" = "Mixing weight (pi):       {round(x$pi, 4)}",
    "*" = "Log-likelihood:           {round(x$log_lik, 2)}",
    "*" = "Converged: {x$converged} ({x$iterations} iterations, n = {x$n})"
  ))
  invisible(x)
}


# Moment-based initialization ----
.init_prior <- function(rr) {
  med <- stats::median(rr)
  low <- rr[rr <= med]
  high <- rr[rr > med]

  # Fit Gamma to each half
  fit_low <- gamma_shape_mle(low)
  fit_high <- gamma_shape_mle(high)

  list(
    alpha1 = fit_low$alpha,
    beta1 = fit_low$beta,
    alpha2 = fit_high$alpha,
    beta2 = fit_high$beta,
    pi = 0.5
  )
}
