# Reporting Odds Ratio (ROR) ----

#' Reporting Odds Ratio (van Puijenbroek 2002)
#'
#' Computes ROR and its 95% CI using the log-normal approximation. Signal
#' criterion: the lower 95% CI exceeds 1 AND observed count a >= min_count.
#'
#' Formulas:
#' \deqn{ROR = \frac{ad}{bc}}
#' \deqn{SE^2(\log ROR) = 1/a + 1/b + 1/c + 1/d}
#'
#' @param oe_data A tibble from [compute_observed_expected()] or with `drug`,
#'   `event`, `observed` columns.
#' @param min_count Minimum `a` (observed count) for signalling. Default 3.
#'
#' @return The input tibble with added columns: `ror`, `ror_lci`, `ror_uci`
#'   (95% CI), `is_signal_ror` (logical).
#'
#' @references
#' van Puijenbroek EP, Bate A, Leufkens HGM, Lindquist M, Orre R, Egberts ACG
#' (2002). "A comparison of measures of disproportionality for signal detection
#' in spontaneous reporting systems for adverse drug reactions."
#' *Pharmacoepidemiology and Drug Safety* 11(1):3-10.
#'
#' Rothman KJ, Lanes S, Sacks ST (2004). "The reporting odds ratio and its
#' advantages over the proportional reporting ratio."
#' *Pharmacoepidemiology and Drug Safety* 13(8):519-523.
#'
#' @export
compute_ror <- function(oe_data, min_count = 3) {
  if (!all(c("a", "b", "c", "d") %in% names(oe_data))) {
    oe_data <- compute_contingency_2x2(oe_data)
  }

  a <- as.double(oe_data$a)
  b <- as.double(oe_data$b)
  c_ <- as.double(oe_data$c)
  d <- as.double(oe_data$d)

  ror <- (a * d) / (b * c_)
  log_ror_var <- 1 / a + 1 / b + 1 / c_ + 1 / d
  log_ror_se <- sqrt(log_ror_var)
  log_ror <- log(ror)
  ror_lci <- exp(log_ror - 1.96 * log_ror_se)
  ror_uci <- exp(log_ror + 1.96 * log_ror_se)

  is_signal_ror <- !is.na(ror_lci) & a >= min_count & ror_lci > 1

  oe_data |>
    dplyr::mutate(
      ror = ror,
      ror_lci = ror_lci,
      ror_uci = ror_uci,
      is_signal_ror = is_signal_ror
    )
}
