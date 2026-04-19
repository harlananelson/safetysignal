test_that("detect_all_methods runs all four methods and returns expected columns", {
  set.seed(42)
  # Build a realistic-enough dataset so EM prior fit converges
  data <- tibble::tibble(
    drug = sample(LETTERS[1:10], 500, replace = TRUE),
    event = sample(paste0("E", 1:20), 500, replace = TRUE)
  )
  result <- detect_all_methods(data, drug, event, verbose = FALSE)

  # All method outputs present
  expect_true("eb05" %in% names(result))
  expect_true("prr" %in% names(result))
  expect_true("ror" %in% names(result))
  expect_true("ic" %in% names(result))
  # Signal flags
  expect_true("is_signal_gps" %in% names(result))
  expect_true("is_signal_prr" %in% names(result))
  expect_true("is_signal_ror" %in% names(result))
  expect_true("is_signal_ic" %in% names(result))
  # Aggregate
  expect_true("n_methods_flagged" %in% names(result))
  expect_true("is_signal_any" %in% names(result))
})

test_that("detect_all_methods subset of methods works", {
  set.seed(42)
  data <- tibble::tibble(
    drug = sample(LETTERS[1:5], 200, replace = TRUE),
    event = sample(paste0("E", 1:5), 200, replace = TRUE)
  )
  result <- detect_all_methods(data, drug, event,
                               methods = c("prr", "ror"),
                               verbose = FALSE)
  expect_true("prr" %in% names(result))
  expect_true("ror" %in% names(result))
  expect_false("eb05" %in% names(result))
  expect_false("ic" %in% names(result))
})

test_that("detect_all_methods accepts pre-aggregated OE table", {
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(20, 80, 10, 890)
  )
  result <- detect_all_methods(oe,
                               methods = c("prr", "ror", "ic"),
                               verbose = FALSE)
  expect_equal(nrow(result), 4)
  expect_true("prr" %in% names(result))
})

test_that("n_methods_flagged correctly counts flagged methods per pair", {
  set.seed(123)
  # Construct a dataset where strong signals should be flagged by multiple methods
  data <- tibble::tibble(
    drug = c(rep("sig_drug", 50),
             sample(LETTERS[2:10], 450, replace = TRUE)),
    event = c(rep("sig_event", 50),
              sample(paste0("E", 2:10), 450, replace = TRUE))
  )
  result <- detect_all_methods(data, drug, event, verbose = FALSE)
  sig <- result[result$drug == "sig_drug" & result$event == "sig_event", ]

  # This synthetic strong signal should hit PRR, ROR, and IC at minimum
  expect_gte(sig$n_methods_flagged, 2)
  expect_true(sig$is_signal_any)
})
