# Observed / Expected table ----

#' Compute observed and expected counts for drug-event pairs
#'
#' Builds the observed/expected table from contingency data using the
#' independence assumption: `E_ij = (n_i. * n_.j) / n_..`
#'
#' @param data A data frame or tibble of drug-event reports.
#' @param drug_col Column identifying the drug/vaccine (tidy-select).
#' @param event_col Column identifying the adverse event (tidy-select).
#' @param count_col Optional column with pre-aggregated counts. If `NULL`,
#'   each row is treated as one report.
#'
#' @return A tibble with columns: `drug`, `event`, `observed`, `expected`, `rr`.
#'
#' @export
compute_observed_expected <- function(data, drug_col, event_col, count_col = NULL) {
  drug_col <- rlang::enquo(drug_col)
  event_col <- rlang::enquo(event_col)
  count_col <- rlang::enquo(count_col)

  # Aggregate counts
  if (rlang::quo_is_null(count_col)) {
    counts <- data |>
      dplyr::count(
        drug = !!drug_col,
        event = !!event_col,
        name = "observed"
      )
  } else {
    counts <- data |>
      dplyr::group_by(
        drug = !!drug_col,
        event = !!event_col
      ) |>
      dplyr::summarize(observed = sum(!!count_col), .groups = "drop")
  }

  # Marginals
  n_total <- sum(counts$observed)
  drug_totals <- counts |>
    dplyr::group_by(.data$drug) |>
    dplyr::summarize(n_drug = sum(.data$observed), .groups = "drop")
  event_totals <- counts |>
    dplyr::group_by(.data$event) |>
    dplyr::summarize(n_event = sum(.data$observed), .groups = "drop")

  # Expected under independence. Coerce to double before multiplying: at
  # full-FAERS scale (millions of reports) n_drug * n_event overflows int32.
  counts |>
    dplyr::left_join(drug_totals, by = "drug") |>
    dplyr::left_join(event_totals, by = "event") |>
    dplyr::mutate(
      expected = (as.double(.data$n_drug) * as.double(.data$n_event)) / as.double(n_total),
      rr = as.double(.data$observed) / .data$expected
    ) |>
    dplyr::select("drug", "event", "observed", "expected", "rr")
}
