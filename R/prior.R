# Prior fitting for the 2-component Gamma mixture ----

#' Fit the 2-component Gamma mixture prior
#'
#' Estimates the hyperparameters of the GPS prior
#' `pi * Gamma(alpha1, beta1) + (1 - pi) * Gamma(alpha2, beta2)` on the
#' relative reporting rate lambda.
#'
#' **Method `"marginal"` (default for observed/expected tables)** is the
#' empirical Bayes fit of DuMouchel (1999): hyperparameters maximise the
#' marginal likelihood of the observed counts,
#' `N | E ~ pi * NB(alpha1, beta1 / (beta1 + E)) + (1 - pi) * NB(alpha2, beta2 / (beta2 + E))`,
#' zero-truncated by default because contingency tables only hold `N >= 1`.
#' Truncation is applied per component, as in DuMouchel (1999) and
#' `openEBGM::negLLsquash()`, so fitted hyperparameters reproduce openEBGM's
#' on the same data and [compute_posterior()] (which matches `openEBGM::Qn()`)
#' completes the reference pipeline.
#' Because each pair contributes through its count *and* its expected count,
#' pairs with tiny `E` carry little information and cannot drag the prior's
#' tail out; the posterior then shrinks low-count pairs toward the bulk of
#' the data. This is what keeps EB05 bounded for products with few reports.
#'
#' **Method `"rr_em"` (legacy)** fits the mixture directly to the observed
#' ratios `N / E` by EM. It ignores Poisson noise in `N / E`, so in sparse
#' data the ratios of single reports with tiny `E` (often in the thousands)
#' inflate the signal component and EB05 for low-count pairs. Kept only for
#' reproducing old outputs and for numeric-vector input.
#'
#' @param data A tibble from [compute_observed_expected()] (`observed` and
#'   `expected` columns for `"marginal"`; `rr` for `"rr_em"`), or a numeric
#'   vector of ratios (legacy `"rr_em"` only).
#' @param method `"marginal"` or `"rr_em"`. Default: `"marginal"` when `data`
#'   has `observed` and `expected`, otherwise `"rr_em"`.
#' @param truncated For `"marginal"`: use the zero-truncated likelihood.
#'   Default `TRUE` (contingency tables omit zero cells).
#' @param squash For `"marginal"`: collapse pairs with identical `observed`
#'   and near-identical `expected` into weighted points before optimising
#'   (DuMouchel & Pregibon 2001 "data squashing"). Default `TRUE`. Changes the
#'   fit negligibly and makes full-FAERS fits fast.
#' @param squash_bin_width Width of the `log(expected)` bins used by `squash`.
#'   Default `0.02` (about 2% relative difference in `E`).
#' @param tol Convergence tolerance. Default `1e-6`.
#' @param max_iter Maximum iterations. Default 200 (`"rr_em"`) or 1000
#'   optimiser iterations per start (`"marginal"`).
#' @param init Optional named list `list(alpha1, beta1, alpha2, beta2, pi)`
#'   of starting values. For `"marginal"` it is tried in addition to the
#'   default starts.
#'
#' @return A list of class `"gps_prior"` with components `alpha1`, `beta1`,
#'   `alpha2`, `beta2`, `pi`, `log_lik`, `converged`, `iterations`, `n`
#'   (pairs used) and `method`.
#'
#' @references DuMouchel W (1999). Bayesian data mining in large frequency
#'   tables, with an application to the FDA spontaneous reporting system.
#'   *The American Statistician* 53(3):177-190.
#'
#' @export
fit_prior <- function(data, method = NULL, truncated = TRUE, squash = TRUE,
                      squash_bin_width = 0.02, tol = 1e-6, max_iter = NULL,
                      init = NULL) {
  has_counts <- is.data.frame(data) &&
    all(c("observed", "expected") %in% names(data))
  if (is.null(method)) {
    method <- if (has_counts) "marginal" else "rr_em"
  }
  method <- match.arg(method, c("marginal", "rr_em"))

  if (method == "marginal") {
    if (!has_counts) {
      cli::cli_abort(c(
        "Method {.val marginal} needs {.field observed} and {.field expected} columns.",
        "i" = "Pass the output of {.fn compute_observed_expected}."
      ))
    }
    return(.fit_prior_marginal(
      observed = data$observed, expected = data$expected,
      truncated = truncated, squash = squash,
      squash_bin_width = squash_bin_width, tol = tol,
      max_iter = if (is.null(max_iter)) 1000L else max_iter, init = init
    ))
  }

  # ---- legacy: EM on observed ratios ----
  if (is.data.frame(data)) {
    if (!"rr" %in% names(data)) {
      cli::cli_abort("{.arg data} must have a column named {.field rr}.")
    }
    rr <- data$rr
  } else {
    rr <- as.numeric(data)
  }
  .fit_prior_rr_em(rr, tol = tol,
                   max_iter = if (is.null(max_iter)) 200L else max_iter,
                   init = init)
}

.fit_prior_rr_em <- function(rr, tol, max_iter, init) {
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
  result$method <- "rr_em"
  class(result) <- "gps_prior"
  result
}


# Empirical Bayes (marginal likelihood) fit ----

# Bounds on the hyperparameters. Without them the zero-truncated likelihood
# can run to degenerate solutions (a component with shape -> 0 behaving as a
# point mass), which fit the counts no better but make posteriors unstable.
.PRIOR_BOUNDS <- list(
  alpha = c(1e-3, 1e3),
  beta  = c(1e-3, 1e3),
  pi    = c(1e-3, 1 - 1e-3)
)

# log(1 - exp(x)) for x <= 0, accurate at both ends (Maechler 2012).
.log1mexp <- function(x) {
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

# Collapse (observed, expected) into weighted points. Pairs share a point when
# they have the same count and log(E) falls in the same narrow bin; E for the
# point is the mean E of its members.
.squash_counts <- function(observed, expected, bin_width) {
  bin <- round(log(expected) / bin_width)
  key <- paste(observed, bin, sep = "_")
  w <- as.numeric(table(key)[unique(key)])
  idx <- match(unique(key), key)
  e_mean <- tapply(expected, key, mean)[unique(key)]
  list(observed = observed[idx], expected = as.numeric(e_mean), weight = w)
}

# Negative marginal log-likelihood on the unconstrained scale
# theta = (log a1, log b1, log a2, log b2, logit pi).
#
# Truncation follows DuMouchel (1999) and openEBGM: each component's negative
# binomial is truncated to N >= 1 separately, then mixed with weight pi,
#   f*(n) = pi * f1(n) / (1 - f1(0)) + (1 - pi) * f2(n) / (1 - f2(0)).
# Posterior weights downstream (compute_posterior) use the untruncated f_k with
# this pi, exactly as openEBGM::Qn does, so the full chain reproduces the
# reference implementation.
.marginal_nll <- function(theta, n, e, w, truncated) {
  a1 <- exp(theta[1]); b1 <- exp(theta[2])
  a2 <- exp(theta[3]); b2 <- exp(theta[4])
  p  <- stats::plogis(theta[5])
  q1 <- b1 / (b1 + e); q2 <- b2 / (b2 + e)
  l1 <- stats::dnbinom(n, size = a1, prob = q1, log = TRUE)
  l2 <- stats::dnbinom(n, size = a2, prob = q2, log = TRUE)
  if (truncated) {
    # log P(N >= 1) = log(1 - q^alpha)
    l1 <- l1 - .log1mexp(pmin(a1 * log(q1), -1e-300))
    l2 <- l2 - .log1mexp(pmin(a2 * log(q2), -1e-300))
  }
  l1 <- l1 + log(p); l2 <- l2 + log1p(-p)
  m <- pmax(l1, l2)
  ll <- m + log(exp(l1 - m) + exp(l2 - m))
  out <- -sum(w * ll)
  if (!is.finite(out)) out <- .Machine$double.xmax / 1e10
  out
}

.fit_prior_marginal <- function(observed, expected, truncated, squash,
                                squash_bin_width, tol, max_iter, init) {
  ok <- is.finite(observed) & is.finite(expected) & expected > 0 &
    (if (truncated) observed >= 1 else observed >= 0)
  observed <- observed[ok]
  expected <- expected[ok]
  n_pairs <- length(observed)
  if (n_pairs < 10) {
    cli::cli_abort("Need at least 10 pairs with positive expected counts to fit the prior, got {n_pairs}.")
  }

  pts <- if (squash) {
    .squash_counts(observed, expected, squash_bin_width)
  } else {
    list(observed = observed, expected = expected, weight = rep(1, n_pairs))
  }

  # Starting values: DuMouchel's (1999) published FAERS starts plus a few
  # spread-out alternatives; best converged optimum wins.
  starts <- list(
    list(alpha1 = 0.2, beta1 = 0.1, alpha2 = 2, beta2 = 4, pi = 1 / 3),
    list(alpha1 = 2, beta1 = 2, alpha2 = 1, beta2 = 0.2, pi = 0.8),
    list(alpha1 = 5, beta1 = 5, alpha2 = 0.5, beta2 = 0.05, pi = 0.95),
    list(alpha1 = 1, beta1 = 1, alpha2 = 3, beta2 = 1, pi = 0.5)
  )
  if (!is.null(init)) starts <- c(list(init), starts)

  to_theta <- function(s) {
    c(log(s$alpha1), log(s$beta1), log(s$alpha2), log(s$beta2),
      stats::qlogis(s$pi))
  }
  lower <- c(log(.PRIOR_BOUNDS$alpha[1]), log(.PRIOR_BOUNDS$beta[1]),
             log(.PRIOR_BOUNDS$alpha[1]), log(.PRIOR_BOUNDS$beta[1]),
             stats::qlogis(.PRIOR_BOUNDS$pi[1]))
  upper <- c(log(.PRIOR_BOUNDS$alpha[2]), log(.PRIOR_BOUNDS$beta[2]),
             log(.PRIOR_BOUNDS$alpha[2]), log(.PRIOR_BOUNDS$beta[2]),
             stats::qlogis(.PRIOR_BOUNDS$pi[2]))

  best <- NULL
  for (s in starts) {
    th0 <- pmin(pmax(to_theta(s), lower), upper)
    # nlminb (PORT) with box constraints: the optimiser openEBGM uses, and
    # markedly more reliable than L-BFGS-B on this flat, correlated surface.
    fit <- tryCatch(
      stats::nlminb(
        th0, .marginal_nll,
        n = pts$observed, e = pts$expected, w = pts$weight,
        truncated = truncated,
        lower = lower, upper = upper,
        control = list(iter.max = max_iter, eval.max = 2L * max_iter,
                       rel.tol = tol * 1e-4)
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || !is.finite(fit$objective)) next
    if (is.null(best) || fit$objective < best$objective) best <- fit
  }
  if (is.null(best)) {
    cli::cli_abort("Marginal-likelihood prior fit failed from every starting value.")
  }

  th <- best$par
  a1 <- exp(th[1]); b1 <- exp(th[2]); a2 <- exp(th[3]); b2 <- exp(th[4])
  p <- stats::plogis(th[5])

  # Component 1 = background (smaller mean)
  result <- if (a1 / b1 > a2 / b2) {
    list(alpha1 = a2, beta1 = b2, alpha2 = a1, beta2 = b1, pi = 1 - p)
  } else {
    list(alpha1 = a1, beta1 = b1, alpha2 = a2, beta2 = b2, pi = p)
  }
  at_bound <- any(abs(th - lower) < 1e-6 | abs(th - upper) < 1e-6)
  if (at_bound) {
    cli::cli_warn(c(
      "Prior fit reached a hyperparameter bound.",
      "i" = "Check {.fn check_prior_calibration}; a bound often means one mixture component is unnecessary."
    ))
  }

  result$log_lik <- -best$objective
  result$converged <- best$convergence == 0
  result$iterations <- unname(best$iterations)
  result$n <- n_pairs
  result$n_points <- length(pts$observed)
  result$truncated <- truncated
  result$at_bound <- at_bound
  result$method <- "marginal"
  class(result) <- "gps_prior"
  result
}


#' Check calibration of a fitted GPS prior
#'
#' Empirical Bayes model check: compares how many pairs have each observed
#' count `k` with the number the fitted prior predicts,
#' `sum_ij P(N_ij = k | E_ij, prior)` (zero-truncated if the prior was), and
#' reports how often pairs with very few reports reach the signal threshold.
#' A well-calibrated prior predicts the count distribution closely, and pairs
#' with one report essentially never reach EB05 >= 2.
#'
#' @param oe_data A tibble from [compute_observed_expected()].
#' @param prior A `gps_prior` from [fit_prior()].
#' @param max_count Largest count tabulated separately; larger counts are
#'   pooled. Default 5.
#' @param threshold EB05 signal threshold for the low-count check. Default 2.
#' @param low_count Pairs with `observed <= low_count` are checked against
#'   `threshold`. Default 1.
#'
#' @return A list with `count_fit` (tibble of `observed_count`, `n_pairs`,
#'   `predicted_pairs`, `ratio`) and `low_count` (tibble with the number of
#'   low-count pairs, how many have EB05 >= threshold, and the maximum EB05).
#'
#' @export
check_prior_calibration <- function(oe_data, prior, max_count = 5,
                                    threshold = 2, low_count = 1) {
  if (!inherits(prior, "gps_prior")) {
    cli::cli_abort("{.arg prior} must be a {.cls gps_prior} object from {.fn fit_prior}.")
  }
  n <- oe_data$observed
  e <- oe_data$expected
  truncated <- isTRUE(prior$truncated) || prior$method == "rr_em"

  # Per-component predictive P(N = k), truncated per component when the
  # prior was fit on N >= 1 data (same likelihood as the fit).
  comp <- function(k, a, b) {
    d <- stats::dnbinom(k, a, b / (b + e))
    if (truncated) d / (1 - stats::dnbinom(0, a, b / (b + e))) else d
  }
  pk <- function(k) {
    prior$pi * comp(k, prior$alpha1, prior$beta1) +
      (1 - prior$pi) * comp(k, prior$alpha2, prior$beta2)
  }
  ks <- seq_len(max_count)
  p_k <- lapply(ks, pk)
  pred <- vapply(p_k, sum, numeric(1))
  obs <- vapply(ks, function(k) sum(n == k), numeric(1))
  p_tail <- 1 - Reduce(`+`, p_k) - if (truncated) 0 else pk(0)
  obs_all <- c(obs, sum(n > max_count))
  pred_all <- c(pred, sum(pmax(p_tail, 0)))
  count_fit <- tibble::tibble(
    observed_count = c(as.character(ks), paste0(">", max_count)),
    n_pairs = obs_all,
    predicted_pairs = pred_all,
    ratio = ifelse(pred_all > 0, obs_all / pred_all, NA_real_)
  )

  low <- oe_data[oe_data$observed <= low_count, , drop = FALSE]
  if (nrow(low) > 0) {
    low <- low |>
      compute_posterior(prior) |>
      posterior_percentile(percentile = 0.05, col_name = "eb05")
    low_summary <- tibble::tibble(
      n_pairs = nrow(low),
      n_at_or_above_threshold = sum(low$eb05 >= threshold, na.rm = TRUE),
      max_eb05 = max(low$eb05, na.rm = TRUE)
    )
  } else {
    low_summary <- tibble::tibble(n_pairs = 0L, n_at_or_above_threshold = 0L,
                                  max_eb05 = NA_real_)
  }
  list(count_fit = count_fit, low_count = low_summary)
}


#' @export
print.gps_prior <- function(x, ...) {
  cli::cli_h3("GPS 2-Component Gamma Mixture Prior")
  cli::cli_bullets(c(
    "*" = "Component 1 (background): Gamma({round(x$alpha1, 4)}, {round(x$beta1, 4)}) — mean = {round(x$alpha1 / x$beta1, 4)}",
    "*" = "Component 2 (signal):     Gamma({round(x$alpha2, 4)}, {round(x$beta2, 4)}) — mean = {round(x$alpha2 / x$beta2, 4)}",
    "*" = "Mixing weight (pi):       {round(x$pi, 4)}",
    "*" = "Log-likelihood:           {round(x$log_lik, 2)}",
    "*" = "Fit method:               {if (is.null(x$method)) 'rr_em' else x$method}",
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
