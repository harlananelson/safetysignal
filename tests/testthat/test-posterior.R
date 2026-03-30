test_that("compute_posterior adds expected columns", {
  set.seed(42)
  oe <- tibble::tibble(
    drug = c("A", "B", "C"),
    event = c("X", "X", "Y"),
    observed = c(10, 2, 50),
    expected = c(5, 8, 10),
    rr = c(2, 0.25, 5)
  )
  prior <- structure(
    list(alpha1 = 2, beta1 = 4, alpha2 = 1, beta2 = 0.5,
         pi = 0.7, log_lik = -100, converged = TRUE, iterations = 10, n = 100),
    class = "gps_prior"
  )

  post <- compute_posterior(oe, prior)

  expect_true(all(c("alpha1_post", "beta1_post", "alpha2_post", "beta2_post", "q_post") %in% names(post)))
  expect_equal(nrow(post), 3)

  # Posterior shape = prior shape + observed
  expect_equal(post$alpha1_post, c(12, 4, 52))
  expect_equal(post$alpha2_post, c(11, 3, 51))

  # Posterior rate = prior rate + expected
  expect_equal(post$beta1_post, c(9, 12, 14))
  expect_equal(post$beta2_post, c(5.5, 8.5, 10.5))

  # Posterior weights between 0 and 1

  expect_true(all(post$q_post >= 0 & post$q_post <= 1))
})

test_that("compute_posterior rejects non-gps_prior", {
  oe <- tibble::tibble(drug = "A", event = "X", observed = 5, expected = 3, rr = 5/3)
  expect_error(compute_posterior(oe, list(alpha1 = 1)), "gps_prior")
})
