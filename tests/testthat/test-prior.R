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
