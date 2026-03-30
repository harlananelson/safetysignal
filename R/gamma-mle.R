# Gamma shape MLE via Newton-Raphson ----

#' Newton-Raphson MLE for Gamma shape parameter
#'
#' Estimates the shape parameter of a Gamma distribution via Newton-Raphson
#' on the score equation: `log(alpha) - digamma(alpha) = log(mean(x)) - mean(log(x))`.
#'
#' @param x Numeric vector of positive observations.
#' @param weights Optional numeric vector of weights (for EM M-step).
#' @param tol Convergence tolerance. Default `1e-8`.
#' @param max_iter Maximum iterations. Default 100.
#'
#' @return A list with components:
#'   - `alpha`: Estimated shape parameter.
#'   - `beta`: Estimated rate parameter (`alpha / weighted_mean`).
#'   - `converged`: Logical.
#'   - `iterations`: Number of iterations used.
#'
#' @export
gamma_shape_mle <- function(x, weights = NULL, tol = 1e-8, max_iter = 100) {
  if (length(x) == 0) {
    cli::cli_abort("{.arg x} must have at least one element.")
  }
  if (any(x <= 0)) {
    cli::cli_abort("{.arg x} must contain only positive values.")
  }

  if (is.null(weights)) {
    weights <- rep(1, length(x))
  }
  weights <- weights / sum(weights)

  # Sufficient statistics
  log_x <- log(x)
  mean_x <- sum(weights * x)
  mean_log_x <- sum(weights * log_x)
  s <- log(mean_x) - mean_log_x

  if (s <= 0) {
    cli::cli_abort("Sufficient statistic s = log(mean(x)) - mean(log(x)) must be positive.")
  }

  # Choi-Wette initial estimate
  alpha <- 0.5 / s

  converged <- FALSE
  for (i in seq_len(max_iter)) {
    # Score: log(alpha) - digamma(alpha) - s
    score <- log(alpha) - digamma(alpha) - s
    # Hessian: 1/alpha - trigamma(alpha)
    hessian <- 1 / alpha - trigamma(alpha)

    step <- score / hessian
    alpha_new <- alpha - step

    # Guard against negative values
    if (alpha_new <= 0) {
      alpha_new <- alpha / 2
    }

    if (abs(alpha_new - alpha) < tol) {
      alpha <- alpha_new
      converged <- TRUE
      break
    }
    alpha <- alpha_new
  }

  beta <- alpha / mean_x

  list(
    alpha = alpha,
    beta = beta,
    converged = converged,
    iterations = i
  )
}
