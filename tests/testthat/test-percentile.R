test_that("posterior_percentile matches qgamma for single component", {
  # When q_post = 1, the mixture is a single Gamma
  post <- tibble::tibble(
    drug = "A", event = "X", observed = 10, expected = 5, rr = 2,
    alpha1_post = 12, beta1_post = 6,
    alpha2_post = 11, beta2_post = 5.5,
    q_post = 1.0  # all weight on component 1
  )

  result <- posterior_percentile(post, percentile = 0.05)
  expected_val <- qgamma(0.05, shape = 12, rate = 6)

  expect_equal(result$eb05, expected_val, tolerance = 1e-6)
})

test_that("posterior_percentile ordering: EB05 < EB95", {
  post <- tibble::tibble(
    drug = c("A", "B"), event = c("X", "Y"),
    observed = c(10, 20), expected = c(5, 3),
    rr = c(2, 6.67),
    alpha1_post = c(12, 22), beta1_post = c(9, 7),
    alpha2_post = c(11, 21), beta2_post = c(5.5, 3.5),
    q_post = c(0.6, 0.4)
  )

  result <- post |>
    posterior_percentile(0.05) |>
    posterior_percentile(0.95)

  expect_true(all(result$eb05 < result$eb95))
})

test_that("parallel path is method-preserving (identical to serial)", {
  skip_on_os("windows")
  old <- options(safetysignal.cores = 1L)
  on.exit(options(old), add = TRUE)

  set.seed(42)
  n <- 3000L # above .SS_PARALLEL_MIN, so the mclapply path is exercised
  post <- tibble::tibble(
    q_post = runif(n, 0.3, 0.9),
    alpha1_post = runif(n, 1, 50), beta1_post = runif(n, 1, 30),
    alpha2_post = runif(n, 1, 20), beta2_post = runif(n, 1, 20)
  )

  serial <- posterior_percentile(post, 0.05)$eb05
  options(safetysignal.cores = 2L)
  parallel <- posterior_percentile(post, 0.05)$eb05

  expect_equal(serial, parallel)
})

test_that("posterior_percentile rejects invalid percentile", {
  post <- tibble::tibble(
    drug = "A", event = "X", observed = 1, expected = 1, rr = 1,
    alpha1_post = 2, beta1_post = 2, alpha2_post = 1, beta2_post = 1,
    q_post = 0.5
  )
  expect_error(posterior_percentile(post, percentile = 0), "between 0 and 1")
  expect_error(posterior_percentile(post, percentile = 1), "between 0 and 1")
})
