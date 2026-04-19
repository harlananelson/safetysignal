# Proportional Reporting Ratio (PRR) ----

#' Proportional Reporting Ratio (Evans 2001)
#'
#' Computes PRR, its 95% CI (log-normal approximation), and Yates-corrected
#' chi-squared for each drug-event pair. The standard MHRA signalling
#' criterion is PRR >= 2 AND chi-square >= 4 AND a >= 3.
#'
#' Formulas:
#' \deqn{PRR = \frac{a/(a+b)}{c/(c+d)}}
#' \deqn{SE^2(\log PRR) = 1/a - 1/(a+b) + 1/c - 1/(c+d)}
#' \deqn{\chi^2_{Yates} = \frac{(|ad-bc| - N/2)^2 N}{(a+b)(c+d)(a+c)(b+d)}}
#'
#' @param oe_data A tibble from [compute_observed_expected()] or with `drug`,
#'   `event`, `observed` columns. The 2x2 cells will be computed via
#'   [compute_contingency_2x2()] if not already present.
#' @param min_count Minimum `a` (observed count) for signalling. Default 3.
#' @param prr_threshold PRR threshold for signalling. Default 2.
#' @param chisq_threshold Yates chi-square threshold for signalling. Default 4.
#'
#' @return A tibble with the input plus: `prr`, `prr_lci`, `prr_uci` (95% CI),
#'   `prr_chisq` (Yates), `is_signal_prr` (logical).
#'
#' @references
#' Evans SJW, Waller PC, Davis S (2001). "Use of proportional reporting
#' ratios (PRRs) for signal generation from spontaneous adverse drug
#' reaction reports." *Pharmacoepidemiology and Drug Safety* 10(6):483-486.
#'
#' @export
compute_prr <- function(oe_data,
                        min_count = 3,
                        prr_threshold = 2,
                        chisq_threshold = 4) {
  if (!all(c("a", "b", "c", "d") %in% names(oe_data))) {
    oe_data <- compute_contingency_2x2(oe_data)
  }

  # Coerce to double to avoid integer overflow in the Yates chi-square
  # denominator (a+b)*(c+d)*(a+c)*(b+d), which can exceed .Machine$integer.max
  a <- as.double(oe_data$a)
  b <- as.double(oe_data$b)
  c_ <- as.double(oe_data$c)
  d <- as.double(oe_data$d)
  n_total <- a + b + c_ + d

  # PRR
  prr_num <- a / (a + b)
  prr_den <- c_ / (c_ + d)
  prr <- prr_num / prr_den

  # Log-normal CI
  log_prr_var <- 1 / a - 1 / (a + b) + 1 / c_ - 1 / (c_ + d)
  log_prr_se <- sqrt(log_prr_var)
  log_prr <- log(prr)
  prr_lci <- exp(log_prr - 1.96 * log_prr_se)
  prr_uci <- exp(log_prr + 1.96 * log_prr_se)

  # Yates-corrected chi-square
  chisq_num <- (abs(a * d - b * c_) - n_total / 2)^2 * n_total
  chisq_den <- (a + b) * (c_ + d) * (a + c_) * (b + d)
  prr_chisq <- chisq_num / chisq_den

  # Signal flag (MHRA criterion)
  is_signal_prr <- !is.na(prr) &
    !is.na(prr_chisq) &
    a >= min_count &
    prr >= prr_threshold &
    prr_chisq >= chisq_threshold

  oe_data |>
    dplyr::mutate(
      prr = prr,
      prr_lci = prr_lci,
      prr_uci = prr_uci,
      prr_chisq = prr_chisq,
      is_signal_prr = is_signal_prr
    )
}
