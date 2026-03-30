test_that("detect_signals flags correctly", {
  data <- tibble::tibble(
    drug = c("A", "B", "C", "D"),
    event = rep("X", 4),
    eb05 = c(0.5, 1.2, 3.0, 8.5)
  )

  result <- detect_signals(data, threshold = 1)

  expect_equal(result$is_signal, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(
    as.character(result$signal_strength),
    c("none", "weak", "moderate", "strong")
  )
})

test_that("detect_signals errors on missing column", {
  data <- tibble::tibble(drug = "A", event = "X", rr = 2)
  expect_error(detect_signals(data), "eb05")
})

test_that("detect_signals respects custom threshold", {
  data <- tibble::tibble(
    drug = c("A", "B"),
    event = c("X", "X"),
    eb05 = c(1.5, 2.5)
  )

  result <- detect_signals(data, threshold = 2)
  expect_equal(result$is_signal, c(FALSE, TRUE))
})
