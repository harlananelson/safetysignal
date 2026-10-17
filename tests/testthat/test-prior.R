test_that("fit_prior converges on simulated mixture data", {
  set.seed(42)
  # Simulate from known 2-component mixture
  n <- 2000
  z <- rbinom(n, 1, 0.3)  # 30% signal
  rr <- ifelse(z == 0, rgamma(n, shape = 2, rate = 4), rgamma(n, shape = 1, rate = 0.5))
  rr <- rr[rr > 0]

  prior <- fit_prior(rr)

  expect_s3_class(prior, "gps_prior")
  expect_true(prior$converged)
  expect_gt(prior$pi, 0)
  expect_lt(prior$pi, 1)

  # Component 1 (background) should have smaller mean than component 2 (signal)
  expect_lt(prior$alpha1 / prior$beta1, prior$alpha2 / prior$beta2)
})

test_that("fit_prior rejects too few observations", {
  expect_error(fit_prior(c(1, 2, 3)), "at least 10")
})

test_that("fit_prior accepts tibble input", {
  set.seed(42)
  data <- tibble::tibble(rr = rgamma(100, shape = 2, rate = 2))
  prior <- fit_prior(data)
  expect_s3_class(prior, "gps_prior")
})

test_that("fit_prior log-likelihood is finite", {
  set.seed(42)
  rr <- rgamma(500, shape = 1.5, rate = 1.5)
  prior <- fit_prior(rr)
  expect_true(is.finite(prior$log_lik))
})

test_that("print.gps_prior runs without error", {
  set.seed(42)
  prior <- fit_prior(rgamma(100, 2, 2))
  # cli output goes to message, not stdout
  expect_message(print(prior), "GPS")
})

# ---- Empirical Bayes (marginal likelihood) fit ----

simulate_counts <- function(n, pi, a1, b1, a2, b2, e_range = c(0.05, 20),
                            seed = 1) {
  set.seed(seed)
  e <- exp(runif(n, log(e_range[1]), log(e_range[2])))
  z <- runif(n) < pi
  lambda <- ifelse(z, rgamma(n, a1, b1), rgamma(n, a2, b2))
  obs <- rpois(n, lambda * e)
  keep <- obs >= 1
  tibble::tibble(drug = seq_len(n)[keep], event = 1L,
                 observed = obs[keep], expected = e[keep],
                 rr = obs[keep] / e[keep])
}

test_that("fit_prior defaults to the marginal method for observed/expected tables", {
  oe <- simulate_counts(3000, 0.8, 20, 20, 2, 0.5)
  prior <- suppressWarnings(fit_prior(oe))
  expect_s3_class(prior, "gps_prior")
  expect_identical(prior$method, "marginal")
  expect_identical(fit_prior(rgamma(100, 2, 2))$method, "rr_em")
})

test_that("marginal fit recovers known hyperparameters from simulated counts", {
  # Truth: 80% tight null around lambda = 1, 20% signal with mean 4.
  oe <- simulate_counts(20000, 0.8, 20, 20, 2, 0.5, seed = 7)
  prior <- suppressWarnings(fit_prior(oe))
  expect_true(prior$converged)
  expect_equal(prior$alpha1 / prior$beta1, 1, tolerance = 0.1)
  expect_equal(prior$alpha2 / prior$beta2, 4, tolerance = 0.25)
  expect_equal(prior$pi, 0.8, tolerance = 0.1)
})

test_that("marginal prior keeps single-report pairs below EB05 2 in sparse data", {
  # Sparse null: every pair has lambda = 1, many have tiny E. Under the legacy
  # ratio-EM fit, ratios of single reports (often > 1000) inflate the signal
  # component and single-report pairs reach EB05 in the hundreds.
  set.seed(11)
  n_drug <- 3000; n_event <- 1500; n_rep <- 30000
  pd <- (1 / runif(n_drug))^(1 / 0.9); pd <- pd / sum(pd)
  pe <- (1 / runif(n_event))^(1 / 0.9); pe <- pe / sum(pe)
  reports <- tibble::tibble(
    drug = sample.int(n_drug, n_rep, TRUE, pd),
    event = sample.int(n_event, n_rep, TRUE, pe)
  )
  oe <- compute_observed_expected(reports, drug, event)
  expect_gt(mean(oe$observed == 1), 0.5)

  prior <- suppressWarnings(fit_prior(oe))
  ones <- oe[oe$observed == 1, ][seq_len(min(1500, sum(oe$observed == 1))), ]
  eb_new <- ones |> compute_posterior(prior) |>
    posterior_percentile(0.05, "eb05")
  expect_true(all(eb_new$eb05 < 2))

  legacy <- fit_prior(oe, method = "rr_em")
  eb_old <- ones |> compute_posterior(legacy) |>
    posterior_percentile(0.05, "eb05")
  expect_gt(max(eb_old$eb05), 2)  # documents the defect the new default fixes
})

test_that("check_prior_calibration reproduces the count distribution for a correct prior", {
  oe <- simulate_counts(20000, 0.8, 20, 20, 2, 0.5, seed = 3)
  prior <- suppressWarnings(fit_prior(oe))
  cal <- check_prior_calibration(oe, prior, max_count = 3)
  expect_named(cal, c("count_fit", "low_count"))
  expect_equal(sum(cal$count_fit$n_pairs), nrow(oe))
  expect_equal(sum(cal$count_fit$predicted_pairs), nrow(oe), tolerance = 1e-6)
  expect_true(all(abs(cal$count_fit$ratio - 1) < 0.1))
  expect_equal(cal$low_count$n_at_or_above_threshold, 0L)
})

test_that("marginal method needs observed and expected", {
  expect_error(fit_prior(tibble::tibble(rr = rgamma(50, 2, 2)), method = "marginal"),
               "observed")
})

test_that("marginal fit reproduces openEBGM hyperparameters", {
  skip_if_not_installed("openEBGM")
  oe <- simulate_counts(20000, 0.8, 20, 20, 2, 0.5, seed = 7)
  ours <- suppressWarnings(fit_prior(oe, squash = FALSE))
  theta_ours <- c(ours$alpha1, ours$beta1, ours$alpha2, ours$beta2, ours$pi)
  ref <- openEBGM::autoHyper(
    data.frame(N = oe$observed, E = oe$expected),
    # autoHyper needs several starts to agree before it reports convergence
    theta_init = data.frame(alpha1 = c(0.2, 2, 5, 1), beta1 = c(0.1, 2, 5, 1),
                            alpha2 = c(2, 1, 0.5, 3), beta2 = c(4, 0.2, 0.05, 1),
                            p = c(1 / 3, 0.8, 0.95, 0.5)),
    squashed = FALSE, zeroes = FALSE, N_star = 1, conf_ints = FALSE
  )
  expect_equal(theta_ours, unname(ref$estimates), tolerance = 1e-3)
})
