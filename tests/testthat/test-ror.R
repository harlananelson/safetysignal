test_that("compute_ror matches ad/bc on a textbook 2x2", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )

  result <- compute_ror(oe)
  ax <- result[result$drug == "A" & result$event == "X", ]

  # ROR = (a*d)/(b*c) = (20*890)/(80*10) = 22.25
  expected_ror <- (20 * 890) / (80 * 10)
  expect_equal(ax$ror, expected_ror, tolerance = 1e-8)
  expect_gt(ax$ror_lci, 1)
  expect_lt(ax$ror_lci, ax$ror)
  expect_gt(ax$ror_uci, ax$ror)
})

test_that("ROR CI uses log-normal approximation correctly", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  result <- compute_ror(oe)
  ax <- result[result$drug == "A" & result$event == "X", ]

  a <- 20; b <- 80; c_ <- 10; d <- 890
  se <- sqrt(1 / a + 1 / b + 1 / c_ + 1 / d)
  expected_lci <- exp(log(ax$ror) - 1.96 * se)
  expected_uci <- exp(log(ax$ror) + 1.96 * se)
  expect_equal(ax$ror_lci, expected_lci, tolerance = 1e-8)
  expect_equal(ax$ror_uci, expected_uci, tolerance = 1e-8)
})

test_that("ROR signal requires lower CI > 1 AND a >= min_count", {
  oe <- tibble::tibble(
    drug = c("big", "small"),
    event = c("E", "E"),
    observed = c(20, 2)
  )
  # Pad with a comparator drug so margins aren't degenerate
  oe <- tibble::tibble(
    drug = c("big", "big", "small", "small", "comp", "comp"),
    event = c("E", "other", "E", "other", "E", "other"),
    observed = c(20, 80, 2, 50, 10, 1000)
  )
  result <- compute_ror(oe)

  big <- result[result$drug == "big" & result$event == "E", ]
  small <- result[result$drug == "small" & result$event == "E", ]

  expect_true(big$is_signal_ror)
  expect_false(small$is_signal_ror)  # a=2 < min_count=3
})
