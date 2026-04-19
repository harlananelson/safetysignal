# 2x2 contingency table ----

#' Build 2x2 contingency cells from an observed/expected table
#'
#' Given the output of [compute_observed_expected()], reconstruct the
#' four cells of the drug-by-event 2x2 table needed by frequentist
#' disproportionality methods (PRR, ROR) and BCPNN/IC:
#'
#' ```
#'               Event+   Event-   Total
#'   Drug+         a        b      a+b
#'   Drug-         c        d      c+d
#'   Total       a+c      b+d       N
#' ```
#'
#' @param oe_data A tibble from [compute_observed_expected()] with columns
#'   `drug`, `event`, `observed` at minimum. Margins are re-derived from
#'   `observed` to match the 2x2 convention (so the helper is self-contained
#'   and doesn't depend on upstream margin columns).
#'
#' @return The input tibble with added columns `a`, `b`, `c`, `d`, `n_total`.
#'
#' @export
compute_contingency_2x2 <- function(oe_data) {
  required <- c("drug", "event", "observed")
  missing_cols <- setdiff(required, names(oe_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg oe_data} is missing required column{?s}: {.field {missing_cols}}",
      call = rlang::caller_env()
    )
  }

  n_total <- sum(oe_data$observed)

  n_drug <- oe_data |>
    dplyr::group_by(.data$drug) |>
    dplyr::summarize(n_drug = sum(.data$observed), .groups = "drop")

  n_event <- oe_data |>
    dplyr::group_by(.data$event) |>
    dplyr::summarize(n_event = sum(.data$observed), .groups = "drop")

  oe_data |>
    dplyr::left_join(n_drug, by = "drug") |>
    dplyr::left_join(n_event, by = "event") |>
    dplyr::mutate(
      a = .data$observed,
      b = .data$n_drug - .data$a,
      c = .data$n_event - .data$a,
      d = n_total - .data$a - .data$b - .data$c,
      n_total = n_total
    ) |>
    dplyr::select(-"n_drug", -"n_event")
}
