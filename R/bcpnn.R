# BCPNN / Information Component (IC) ----

#' Information Component from Bayesian Confidence Propagation Neural Network
#'
#' Computes the Information Component (IC) and its lower 95% credible bound
#' (IC025) per Bate et al. 1998 / Noren approximation. Signal criterion:
#' IC025 > 0 (i.e., the 95% posterior credibility excludes the null).
#'
#' Formulas (Bate/Noren approximation with weakly-informative priors
#' alpha = beta = gamma = 0.5):
#' \deqn{IC = \log_2 \frac{(a + 0.5)/(N + 1)}{[(a+b+0.5)(a+c+0.5)] / (N+1)^2}}
#' \deqn{V(IC) \approx \frac{1}{(\ln 2)^2}\left[
#'   \frac{N - a + \gamma - 1}{(a+1)(N+\gamma+1)} +
#'   \frac{N - (a+b) + \alpha - 1}{(a+b+1)(N+\alpha+1)} +
#'   \frac{N - (a+c) + \beta - 1}{(a+c+1)(N+\beta+1)}
#' \right]}
#' \deqn{IC_{025} = IC - 1.96 \sqrt{V(IC)}}
#'
#' @param oe_data A tibble from [compute_observed_expected()] or with `drug`,
#'   `event`, `observed` columns.
#'
#' @return The input tibble with added columns: `ic` (point estimate),
#'   `ic025` (lower 95% credible bound), `ic975` (upper 95% credible bound),
#'   `is_signal_ic` (logical, `ic025 > 0`).
#'
#' @references
#' Bate A, Lindquist M, Edwards IR, Olsson S, Orre R, Lansner A, De Freitas RM
#' (1998). "A Bayesian neural network method for adverse drug reaction signal
#' generation." *European Journal of Clinical Pharmacology* 54(4):315-321.
#'
#' Noren GN, Bate A, Orre R, Edwards IR (2006). "Extending the methods used to
#' screen the WHO drug safety database towards analysis of complex associations
#' and improved accuracy for rare events." *Statistics in Medicine* 25(21):
#' 3740-3757.
#'
#' @export
compute_ic <- function(oe_data) {
  if (!all(c("a", "b", "c", "d") %in% names(oe_data))) {
    oe_data <- compute_contingency_2x2(oe_data)
  }

  a <- as.double(oe_data$a)
  b <- as.double(oe_data$b)
  c_ <- as.double(oe_data$c)
  d <- as.double(oe_data$d)
  n <- a + b + c_ + d

  # Bayesian IC with weakly-informative priors (alpha = beta = gamma = 0.5)
  alpha <- 0.5
  beta_ <- 0.5
  gamma <- 0.5

  p_ab <- (a + gamma) / (n + 1)
  p_a_marg <- (a + b + alpha) / (n + 1)
  p_b_marg <- (a + c_ + beta_) / (n + 1)

  ic <- log2(p_ab / (p_a_marg * p_b_marg))

  # Variance approximation (Noren 2006)
  ln2_sq <- log(2)^2
  v1 <- (n - a + gamma - 1) / ((a + 1) * (n + gamma + 1))
  v2 <- (n - (a + b) + alpha - 1) / ((a + b + 1) * (n + alpha + 1))
  v3 <- (n - (a + c_) + beta_ - 1) / ((a + c_ + 1) * (n + beta_ + 1))
  var_ic <- (v1 + v2 + v3) / ln2_sq
  # Clamp negatives (can occur with very small N)
  var_ic <- pmax(var_ic, 0)
  se_ic <- sqrt(var_ic)

  ic025 <- ic - 1.96 * se_ic
  ic975 <- ic + 1.96 * se_ic

  is_signal_ic <- !is.na(ic025) & ic025 > 0

  oe_data |>
    dplyr::mutate(
      ic = ic,
      ic025 = ic025,
      ic975 = ic975,
      is_signal_ic = is_signal_ic
    )
}
