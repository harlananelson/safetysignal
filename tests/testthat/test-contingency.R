test_that("compute_contingency_2x2 produces correct cells for a known table", {
  # 2x2: drug A has 10 reports with event X, 5 with Y
  #      drug B has 3 reports with X, 12 with Y
  # Totals: A+X=10, A+Y=5, B+X=3, B+Y=12, N=30
  # For (A, X): a=10, b=5 (A+other=Y), c=3 (other+X=B+X), d=12 (other+other=B+Y)
  oe <- tibble::tibble(
    drug = c("A", "A", "B", "B"),
    event = c("X", "Y", "X", "Y"),
    observed = c(10, 5, 3, 12)
  )

  ct <- compute_contingency_2x2(oe)

  expect_true(all(c("a", "b", "c", "d", "n_total") %in% names(ct)))

  ax <- ct[ct$drug == "A" & ct$event == "X", ]
  expect_equal(ax$a, 10)
  expect_equal(ax$b, 5)
  expect_equal(ax$c, 3)
  expect_equal(ax$d, 12)
  expect_equal(ax$n_total, 30)

  by_ <- ct[ct$drug == "B" & ct$event == "Y", ]
  expect_equal(by_$a, 12)
  expect_equal(by_$b, 3)
  expect_equal(by_$c, 5)
  expect_equal(by_$d, 10)
})

test_that("a + b + c + d always equals n_total for every row", {
  set.seed(42)
  data <- tibble::tibble(
    drug = sample(LETTERS[1:5], 200, replace = TRUE),
    event = sample(paste0("AE", 1:10), 200, replace = TRUE)
  )
  oe <- compute_observed_expected(data, drug, event)
  ct <- compute_contingency_2x2(oe)

  expect_true(all(ct$a + ct$b + ct$c + ct$d == ct$n_total))
})

test_that("compute_contingency_2x2 errors when required columns missing", {
  bad <- tibble::tibble(x = 1, y = 2)
  expect_error(compute_contingency_2x2(bad), "missing required column")
})
