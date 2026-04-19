test_that("compute_ic returns ic, ic025, ic975 columns", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  result <- compute_ic(oe)
  expect_true(all(c("ic", "ic025", "ic975", "is_signal_ic") %in% names(result)))
  # ic025 < ic < ic975 for each row
  expect_true(all(result$ic025 <= result$ic))
  expect_true(all(result$ic <= result$ic975))
})

test_that("IC is positive when disproportionality is high", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  result <- compute_ic(oe)
  ax <- result[result$drug == "A" & result$event == "X", ]
  expect_gt(ax$ic, 0)
  expect_gt(ax$ic025, 0)  # 95% lower bound excludes null -> signal
  expect_true(ax$is_signal_ic)
})

test_that("IC is near zero when drug and event are independent", {
  # Perfect independence: A & B have proportional (X, Y) splits
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(50, 50, 50, 50)
  )
  result <- compute_ic(oe)
  # All four cells should have ic ~ 0
  expect_true(all(abs(result$ic) < 0.1))
})

test_that("IC point estimate matches Bayesian log2 formula", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  result <- compute_ic(oe)
  ax <- result[result$drug == "A" & result$event == "X", ]

  a <- 20; b <- 80; c_ <- 10; d <- 890; n <- 1000
  alpha <- 0.5; beta_ <- 0.5; gamma <- 0.5
  p_ab <- (a + gamma) / (n + 1)
  p_a_marg <- (a + b + alpha) / (n + 1)
  p_b_marg <- (a + c_ + beta_) / (n + 1)
  expected_ic <- log2(p_ab / (p_a_marg * p_b_marg))
  expect_equal(ax$ic, expected_ic, tolerance = 1e-8)
})
