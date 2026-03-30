test_that("compute_observed_expected produces correct values for known table", {
  # Simple 2x2: drug A has 10 reports with event X, 5 with event Y
  #             drug B has 3 reports with event X, 12 with event Y
  data <- tibble::tibble(
    drug = c(rep("A", 15), rep("B", 15)),
    event = c(rep("X", 10), rep("Y", 5), rep("X", 3), rep("Y", 12))
  )

  oe <- compute_observed_expected(data, drug, event)

  expect_equal(nrow(oe), 4)
  expect_true(all(c("drug", "event", "observed", "expected", "rr") %in% names(oe)))

  # Check A-X: observed = 10, expected = 15 * 13 / 30 = 6.5

  ax <- oe[oe$drug == "A" & oe$event == "X", ]
  expect_equal(ax$observed, 10)
  expect_equal(ax$expected, 6.5)
  expect_equal(ax$rr, 10 / 6.5)
})

test_that("compute_observed_expected works with count column", {
  data <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    n = c(10, 5, 3, 12)
  )

  oe <- compute_observed_expected(data, drug, event, n)
  expect_equal(nrow(oe), 4)

  ax <- oe[oe$drug == "A" & oe$event == "X", ]
  expect_equal(ax$observed, 10)
})

test_that("expected values sum to total observed", {
  set.seed(42)
  data <- tibble::tibble(
    drug = sample(LETTERS[1:5], 200, replace = TRUE),
    event = sample(paste0("AE", 1:10), 200, replace = TRUE)
  )

  oe <- compute_observed_expected(data, drug, event)
  expect_equal(sum(oe$observed), sum(oe$expected), tolerance = 1e-10)
})
