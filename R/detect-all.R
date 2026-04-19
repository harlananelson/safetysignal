# Unified multi-method signal detection ----

#' Run multiple disproportionality methods and combine results
#'
#' Runs any combination of GPS/EBGM (DuMouchel 1999), PRR (Evans 2001),
#' ROR (van Puijenbroek 2002), and BCPNN/IC (Bate 1998) on the same
#' drug-event data, left-joining results into a single wide tibble. Intended
#' as the single entry point for applications and batch compute pipelines.
#'
#' The observed/expected table and the 2x2 contingency are computed once and
#' reused across methods. GPS uses the EM-fitted Gamma mixture prior on the
#' same data (or a pre-fitted prior passed via `prior`).
#'
#' @param data Drug-event report data frame, or a tibble with `drug`, `event`,
#'   `observed` columns (already aggregated).
#' @param drug_col Column identifying the drug (tidy-select). Ignored when
#'   `data` is already an observed/expected table.
#' @param event_col Column identifying the event (tidy-select). Ignored when
#'   `data` is already an observed/expected table.
#' @param count_col Optional pre-aggregated count column.
#' @param methods Character vector of methods to run. Any subset of
#'   `c("gps", "prr", "ror", "ic")`. Default all four.
#' @param gps_threshold EB05 signal threshold for GPS. Default 2 (FDA).
#' @param prr_threshold PRR threshold. Default 2.
#' @param prr_chisq_threshold Yates chi-square threshold. Default 4.
#' @param min_count Minimum `a` for PRR/ROR signalling. Default 3.
#' @param prior Optional pre-fitted `gps_prior` object to use for GPS.
#' @param verbose If `TRUE`, emit cli progress messages. Default `TRUE`.
#'
#' @return A tibble with one row per drug-event pair containing columns from
#'   each requested method plus `n_methods_flagged` (integer count of methods
#'   that flagged the pair) and `is_signal_any` (logical).
#'
#' @export
detect_all_methods <- function(data,
                               drug_col = NULL, event_col = NULL,
                               count_col = NULL,
                               methods = c("gps", "prr", "ror", "ic"),
                               gps_threshold = 2,
                               prr_threshold = 2,
                               prr_chisq_threshold = 4,
                               min_count = 3,
                               prior = NULL,
                               verbose = TRUE) {
  methods <- match.arg(methods, several.ok = TRUE)

  # Accept either already-aggregated OE data or raw reports
  already_oe <- all(c("drug", "event", "observed") %in% names(data))

  if (!already_oe) {
    drug_col <- rlang::enquo(drug_col)
    event_col <- rlang::enquo(event_col)
    count_col <- rlang::enquo(count_col)

    if (verbose) cli::cli_inform("Aggregating observed/expected table...")
    if (rlang::quo_is_null(count_col)) {
      oe <- compute_observed_expected(data, !!drug_col, !!event_col)
    } else {
      oe <- compute_observed_expected(data, !!drug_col, !!event_col, !!count_col)
    }
  } else {
    # Ensure expected/rr exist; compute if only drug/event/observed given
    if (!all(c("expected", "rr") %in% names(data))) {
      # Rebuild via compute_observed_expected-like margins
      n_total <- sum(data$observed)
      drug_totals <- data |>
        dplyr::group_by(.data$drug) |>
        dplyr::summarize(n_drug = sum(.data$observed), .groups = "drop")
      event_totals <- data |>
        dplyr::group_by(.data$event) |>
        dplyr::summarize(n_event = sum(.data$observed), .groups = "drop")
      oe <- data |>
        dplyr::left_join(drug_totals, by = "drug") |>
        dplyr::left_join(event_totals, by = "event") |>
        dplyr::mutate(
          expected = (.data$n_drug * .data$n_event) / n_total,
          rr = .data$observed / .data$expected
        ) |>
        dplyr::select("drug", "event", "observed", "expected", "rr")
    } else {
      oe <- data
    }
  }

  # Build 2x2 cells once (needed by PRR, ROR, IC)
  ct <- compute_contingency_2x2(oe)

  result <- oe

  if ("gps" %in% methods) {
    if (verbose) cli::cli_inform("Running GPS/EBGM...")
    if (is.null(prior)) {
      prior <- fit_prior(oe)
    }
    gps_out <- oe |>
      compute_posterior(prior) |>
      posterior_percentile(percentile = 0.05, col_name = "eb05") |>
      posterior_percentile(percentile = 0.5, col_name = "eb50") |>
      posterior_percentile(percentile = 0.95, col_name = "eb95")
    gps_out$is_signal_gps <- !is.na(gps_out$eb05) & gps_out$eb05 >= gps_threshold
    result <- result |>
      dplyr::left_join(
        gps_out |> dplyr::select("drug", "event", "eb05", "eb50", "eb95",
                                 "is_signal_gps"),
        by = c("drug", "event")
      )
  }

  if ("prr" %in% methods) {
    if (verbose) cli::cli_inform("Running PRR (Evans 2001)...")
    prr_out <- compute_prr(ct,
                           min_count = min_count,
                           prr_threshold = prr_threshold,
                           chisq_threshold = prr_chisq_threshold)
    result <- result |>
      dplyr::left_join(
        prr_out |> dplyr::select("drug", "event", "prr", "prr_lci", "prr_uci",
                                 "prr_chisq", "is_signal_prr"),
        by = c("drug", "event")
      )
  }

  if ("ror" %in% methods) {
    if (verbose) cli::cli_inform("Running ROR (van Puijenbroek 2002)...")
    ror_out <- compute_ror(ct, min_count = min_count)
    result <- result |>
      dplyr::left_join(
        ror_out |> dplyr::select("drug", "event", "ror", "ror_lci", "ror_uci",
                                 "is_signal_ror"),
        by = c("drug", "event")
      )
  }

  if ("ic" %in% methods) {
    if (verbose) cli::cli_inform("Running BCPNN/IC (Bate 1998)...")
    ic_out <- compute_ic(ct)
    result <- result |>
      dplyr::left_join(
        ic_out |> dplyr::select("drug", "event", "ic", "ic025", "ic975",
                                "is_signal_ic"),
        by = c("drug", "event")
      )
  }

  # Aggregate signal flags across methods
  signal_cols <- grep("^is_signal_", names(result), value = TRUE)
  if (length(signal_cols) > 0) {
    flag_matrix <- as.matrix(result[, signal_cols, drop = FALSE])
    # Treat NA as FALSE for counting
    flag_matrix[is.na(flag_matrix)] <- FALSE
    result$n_methods_flagged <- rowSums(flag_matrix)
    result$is_signal_any <- result$n_methods_flagged > 0
  }

  result
}
