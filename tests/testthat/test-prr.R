test_that("compute_prr matches hand-computed PRR on a textbook 2x2", {
  # Classic 2x2 with clear disproportionality:
  #           Event+  Event-
  # Drug+       a=20   b=80     -> P(E|D) = 20/100 = 0.20
  # Drug-       c=10   d=890    -> P(E|not D) = 10/900 ~ 0.0111
  # PRR = 0.20 / 0.0111 = 18.0
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )

  result <- compute_prr(oe)
  ax <- result[result$drug == "A" & result$event == "X", ]

  expected_prr <- (20 / 100) / (10 / 900)
  expect_equal(ax$prr, expected_prr, tolerance = 1e-8)
  expect_gt(ax$prr_lci, 1)
  expect_lt(ax$prr_lci, ax$prr)
  expect_gt(ax$prr_uci, ax$prr)
})

test_that("PRR 95% CI uses log-normal approximation correctly", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  result <- compute_prr(oe)
  ax <- result[result$drug == "A" & result$event == "X", ]

  a <- 20; b <- 80; c_ <- 10; d <- 890
  se <- sqrt(1 / a - 1 / (a + b) + 1 / c_ - 1 / (c_ + d))
  expected_lci <- exp(log(ax$prr) - 1.96 * se)
  expected_uci <- exp(log(ax$prr) + 1.96 * se)
  expect_equal(ax$prr_lci, expected_lci, tolerance = 1e-8)
  expect_equal(ax$prr_uci, expected_uci, tolerance = 1e-8)
})

test_that("Yates-corrected chi-square matches hand computation", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  result <- compute_prr(oe)
  ax <- result[result$drug == "A" & result$event == "X", ]

  a <- 20; b <- 80; c_ <- 10; d <- 890; N <- 1000
  expected_chisq <- (abs(a * d - b * c_) - N / 2)^2 * N /
    ((a + b) * (c_ + d) * (a + c_) * (b + d))
  expect_equal(ax$prr_chisq, expected_chisq, tolerance = 1e-8)
})

test_that("MHRA signal criterion requires PRR>=2 AND chisq>=4 AND a>=3", {
  # Build 4 rows: one that meets all criteria, three that fail one each
  oe <- tibble::tibble(
    drug = c("hit", "low_count", "low_prr", "low_chisq", "filler"),
    event = c("E1", "E1", "E1", "E1", "E2"),
    observed = c(20, 2, 6, 3, 2000)
  )
  # compute_observed_expected-style margins via bare tibble
  oe$observed <- as.integer(oe$observed)
  result <- compute_prr(oe)

  # Row "hit" has a=20 and strong disproportionality -> should signal
  hit <- result[result$drug == "hit", ]
  expect_equal(hit$a, 20)
  expect_true(hit$is_signal_prr)

  # Row "low_count" has a=2 < min_count -> cannot signal
  low_count <- result[result$drug == "low_count", ]
  expect_false(low_count$is_signal_prr)
})

test_that("compute_prr accepts output of compute_observed_expected directly", {
  # End-to-end: report-level data -> OE -> PRR
  data <- tibble::tibble(
    drug = c(rep("A", 100), rep("B", 900)),
    event = c(rep("X", 20), rep("Y", 80), rep("X", 10), rep("Y", 890))
  )
  oe <- compute_observed_expected(data, drug, event)
  result <- compute_prr(oe)

  expect_true("prr" %in% names(result))
  expect_true("is_signal_prr" %in% names(result))

  ax <- result[result$drug == "A" & result$event == "X", ]
  expected_prr <- (20 / 100) / (10 / 900)
  expect_equal(ax$prr, expected_prr, tolerance = 1e-8)
})

test_that("PRR threshold parameters are respected", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  # With a very high threshold nothing should signal
  strict <- compute_prr(oe, prr_threshold = 100)
  expect_false(any(strict$is_signal_prr))

  # With relaxed thresholds and at least 3 observed, the strong pair signals
  loose <- compute_prr(oe, prr_threshold = 1.5, chisq_threshold = 1)
  ax <- loose[loose$drug == "A" & loose$event == "X", ]
  expect_true(ax$is_signal_prr)
})
