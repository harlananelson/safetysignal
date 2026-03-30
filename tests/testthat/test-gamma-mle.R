test_that("gamma_shape_mle recovers known parameters", {
  set.seed(42)
  x <- rgamma(1000, shape = 3, rate = 2)
  fit <- gamma_shape_mle(x)

  expect_true(fit$converged)
  expect_equal(fit$alpha, 3, tolerance = 0.2)
  expect_equal(fit$beta, 2, tolerance = 0.2)
})

test_that("gamma_shape_mle works with weights", {
  set.seed(42)
  x <- rgamma(500, shape = 5, rate = 1)
  w <- runif(500)
  fit <- gamma_shape_mle(x, weights = w)

  expect_true(fit$converged)
  expect_gt(fit$alpha, 0)
  expect_gt(fit$beta, 0)
})

test_that("gamma_shape_mle rejects empty input", {
  expect_error(gamma_shape_mle(numeric(0)), "at least one element")
})

test_that("gamma_shape_mle rejects non-positive values", {
  expect_error(gamma_shape_mle(c(1, -1, 2)), "positive values")
})

test_that("gamma_shape_mle converges for extreme shape", {
  set.seed(42)
  # Very small shape (exponential-like)
  x <- rgamma(500, shape = 0.1, rate = 1)
  fit <- gamma_shape_mle(x)
  expect_true(fit$converged)
  expect_equal(fit$alpha, 0.1, tolerance = 0.1)
})
