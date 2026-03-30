# Signal detection ----

#' Flag drug-event pairs as signals
#'
#' Classifies each pair based on whether the posterior lower bound
#' (e.g., EB05) exceeds the threshold.
#'
#' @param posterior_data A tibble from [posterior_percentile()].
#' @param threshold Signal threshold. Default `1` (EB05 > 1 means 95%
#'   posterior probability that the true RR exceeds 1).
#' @param percentile_col Column containing the percentile to compare.
#'   Default `"eb05"`.
#'
#' @return The input tibble with additional columns:
#'   - `is_signal`: Logical.
#'   - `signal_strength`: Factor with levels `"none"`, `"weak"`,
#'     `"moderate"`, `"strong"`.
#'
#' @export
detect_signals <- function(posterior_data, threshold = 1, percentile_col = "eb05") {
  if (!percentile_col %in% names(posterior_data)) {
    cli::cli_abort(
      "Column {.field {percentile_col}} not found. Run {.fn posterior_percentile} first.",
      call = rlang::caller_env()
    )
  }

  eb_vals <- posterior_data[[percentile_col]]

  posterior_data |>
    dplyr::mutate(
      is_signal = eb_vals >= threshold,
      signal_strength = factor(
        dplyr::case_when(
          eb_vals >= threshold * 4 ~ "strong",
          eb_vals >= threshold * 2 ~ "moderate",
          eb_vals >= threshold ~ "weak",
          .default = "none"
        ),
        levels = c("none", "weak", "moderate", "strong")
      )
    )
}
