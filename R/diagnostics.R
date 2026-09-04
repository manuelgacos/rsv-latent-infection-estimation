# Diagnostics
#
# Provides summary, probability, and plotting helpers for the RSV simulation
# and recovery workflow.


#' Summarize healthcare-visit outcomes
#'
#' Converts an `rsv_data` object into a one-row-per-subject data frame
#' containing birth-timing and healthcare-visit summaries. Optionally computes
#' each subject's model-implied probability of a healthcare visit over the
#' evaluation window.
#'
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases used
#'   to define the age window and compute model-implied visit probabilities.
#' @param visit_age_max Positive integer scalar or `NULL`. Upper bound for
#'   healthcare-visit ages included in the evaluation window. If `NULL`, the
#'   value of `days` is used.
#' @param days Positive integer scalar giving the age-window length. It is used
#'   as `visit_age_max` when that argument is `NULL` and as the evaluation
#'   window when `compute_g = TRUE`. The default is `ncol(model$B_day)`.
#' @param compute_g Logical scalar indicating whether to compute each subject's
#'   model-implied healthcare-visit probability \eqn{g_i}. The default is
#'   `FALSE`.
#' @param beta Numeric vector of length J containing the infection-age spline
#'   coefficients. Required when `compute_g = TRUE`.
#' @param eta Numeric vector of length K containing the healthcare-visit spline
#'   coefficients. Required when `compute_g = TRUE`.
#' @param control An `rsv_control` object controlling numerical checks and
#'   clipping when `compute_g = TRUE`.
#' @param subj_pre_list List of subject-level precomputations or `NULL`.
#'   Required when `compute_g = TRUE`, with one element per subject containing
#'   `lambda_i` and `V_i`.
#' @param month_len Positive numeric scalar giving the number of days used to
#'   map day indices to month-like bins. The default is `365 / 12`.
#' @param nbins Integer scalar at least 2 giving the number of quantile bins
#'   used for `g_decile` when `compute_g = TRUE`. The default is 10.
#'
#' @return A data frame with one row per subject containing:
#' \describe{
#'   \item{\code{i}}{Integer row index.}
#'   \item{\code{id}}{Numeric subject identifier, or `NA` if absent.}
#'   \item{\code{birth_index}}{Integer one-based birth index.}
#'   \item{\code{birth_month}}{Integer month-like bin derived from
#'     `birth_index` and capped at 12.}
#'   \item{\code{birth_decile}}{Integer quantile bin from 1 through 10 based
#'     on `birth_index`.}
#'   \item{\code{visit_age}}{Integer healthcare-visit age, or `NA` if no visit
#'     is recorded.}
#'   \item{\code{visit_in_window}}{Logical indicator for a recorded visit
#'     between days 1 and `visit_age_max`.}
#'   \item{\code{I_visit}}{Integer indicator equal to 1 for an in-window visit
#'     and 0 otherwise.}
#'   \item{\code{visit_month}}{Integer month-like bin for an in-window
#'     `visit_age`, or `NA` otherwise.}
#' }
#'
#' If `compute_g = TRUE`, the data frame also contains:
#' \describe{
#'   \item{\code{g_i}}{Numeric model-implied healthcare-visit probability.}
#'   \item{\code{I_minus_g}}{Numeric difference `I_visit - g_i`.}
#'   \item{\code{abs_I_minus_g}}{Nonnegative numeric absolute difference
#'     `abs(I_visit - g_i)`.}
#'   \item{\code{g_decile}}{Integer quantile bin from 1 through `nbins` based
#'     on `g_i`.}
#' }
#'
#' @details
#' A visit is treated as in-window when `visit_age` is nonmissing and lies
#' between 1 and `visit_age_max`, inclusive.
#'
#' When `compute_g = TRUE`, the model-implied healthcare-visit probability is
#' \deqn{
#'   g_i(\beta,\eta) = \eta^\top U_i(\beta),
#' }
#' where
#' \deqn{
#'   U_i(\beta) = \sum_{m=1}^{D} s(m)Q_i(m;\beta),
#' }
#' and \eqn{Q_i(m;\beta)} is the daily first-infection probability mass over
#' the effective evaluation window of length \eqn{D}.
extract_visit_ages <- function(data,
                               model,
                               visit_age_max = NULL,
                               days = ncol(model$B_day),
                               compute_g = FALSE,
                               beta = NULL,
                               eta  = NULL,
                               control = make_rsv_control(check_bounds = TRUE, warn_on_clip = FALSE),
                               subj_pre_list = NULL,
                               month_len = 365 / 12,
                               nbins = 10L) {
  stopifnot(inherits(data, "rsv_data"))
  stopifnot(is.list(data$subjects), length(data$subjects) >= 1L)
  stopifnot(inherits(model, "rsv_model"))
  
  n <- length(data$subjects)
  days <- as.integer(days)
  stopifnot(days >= 1L)
  
  if (is.null(visit_age_max)) visit_age_max <- days
  visit_age_max <- as.integer(visit_age_max)
  stopifnot(visit_age_max >= 1L)
  
  # Map day indices to capped month-like bins.
  day_to_month <- function(day, month_len, max_month = 12L) {
    ifelse(is.na(day), NA_integer_,
           pmin(max_month, pmax(1L, as.integer(ceiling(day / month_len)))))
  }
  
  # Construct quantile bins while handling tied values.
  quantile_bins <- function(x, nbins = 10L) {
    nbins <- as.integer(nbins)
    stopifnot(nbins >= 2L)
    
    out <- rep(NA_integer_, length(x))
    ok <- is.finite(x)
    if (!any(ok)) return(out)
    
    qs <- stats::quantile(x[ok], probs = seq(0, 1, length.out = nbins + 1L),
                          na.rm = TRUE, type = 7)
    
    # Use rank-based bins when tied values collapse quantile breaks.
    if (any(diff(qs) <= 0)) {
      r <- rank(x[ok], ties.method = "average")
      out[ok] <- pmin(nbins, pmax(1L, ceiling(nbins * r / max(r))))
      return(out)
    }
    
    out[ok] <- as.integer(cut(x[ok], breaks = qs, include.lowest = TRUE, labels = FALSE))
    out
  }
  
  # Extract the subject-level fields used in the summary.
  id <- vapply(data$subjects, function(s) if (!is.null(s$id)) s$id else NA_real_, numeric(1))
  birth_index <- vapply(data$subjects, function(s) as.integer(s$birth_index), integer(1))
  visit_age <- vapply(
    data$subjects,
    function(s) if (is.null(s$visit_age) || is.na(s$visit_age)) NA_integer_ else as.integer(s$visit_age),
    integer(1)
  )
  
  visit_in_window <- !is.na(visit_age) & visit_age >= 1L & visit_age <= visit_age_max
  I_visit <- as.integer(visit_in_window)
  
  df <- data.frame(
    i = seq_len(n),
    id = id,
    birth_index = birth_index,
    birth_month = day_to_month(birth_index, month_len = month_len, max_month = 12L),
    birth_decile = quantile_bins(birth_index, nbins = 10L),
    visit_age = visit_age,
    visit_in_window = visit_in_window,
    I_visit = I_visit,
    visit_month = ifelse(visit_in_window,
                         day_to_month(visit_age, month_len = month_len, max_month = 12L),
                         NA_integer_),
    stringsAsFactors = FALSE
  )
  
  # Optionally compute model-implied healthcare-visit probabilities.
  compute_g <- isTRUE(compute_g)
  if (compute_g) {
    if (is.null(subj_pre_list)) {
      stop("extract_visit_ages: compute_g=TRUE requires subj_pre_list (precomputed via rsv_precompute_subject).")
    }
    stopifnot(is.list(subj_pre_list), length(subj_pre_list) == n)
    stopifnot(is.numeric(beta), is.numeric(eta))
    
    # Precompute coefficient-dependent curves once for all subjects.
    glob <- rsv_precompute_global(beta = beta, eta = eta, model = model, control = control)
    w_vec <- glob$w_vec
    c_vec <- glob$c_vec
    
    # Restrict the evaluation window to dimensions supported by all inputs.
    D_global <- min(days, length(w_vec), length(c_vec), ncol(model$S_day), ncol(model$B_day))
    days_eff <- as.integer(D_global)
    
    w_eff <- w_vec[seq_len(days_eff)]
    S_eff <- model$S_day[, seq_len(days_eff), drop = FALSE]
    
    g_i <- rep(NA_real_, n)
    
    for (k in seq_len(n)) {
      pre <- subj_pre_list[[k]]
      lambda_i <- pre$lambda_i
      V_i <- pre$V_i
      
      Fbar <- Fbar_i(beta = beta, V_i = V_i,
                     include_day0 = TRUE,
                     check_bounds = isTRUE(control$check_bounds),
                     warn_on_clip = isTRUE(control$warn_on_clip))
      
      pi <- pi_from_w(lambda_shift_i = lambda_i, w = w_eff,
                      check_bounds = isTRUE(control$check_bounds),
                      warn_on_clip = isTRUE(control$warn_on_clip))
      
      Q <- Q_i(Fbar_i = Fbar, pi_i = pi,
               check_bounds = isTRUE(control$check_bounds),
               warn_on_clip = isTRUE(control$warn_on_clip))
      
      D <- min(days_eff, length(Q), length(pi), length(Fbar) - 1L)
      Q_eff <- Q[seq_len(D)]
      S_sub <- S_eff[, seq_len(D), drop = FALSE]
      
      U <- U_i(Q_i = Q_eff, S_day = S_sub,
               check_bounds = isTRUE(control$check_bounds),
               warn_on_clip = isTRUE(control$warn_on_clip))
      
      g_i[k] <- sum(eta * U)
    }
    
    df$g_i <- g_i
    df$I_minus_g <- df$I_visit - df$g_i
    df$abs_I_minus_g <- abs(df$I_minus_g)
    df$g_decile <- quantile_bins(df$g_i, nbins = nbins)
  }
  
  df
}


#' Compute first-infection probability by age
#'
#' Computes the model-implied probability that a subject experiences their
#' first RSV infection by a specified age.
#'
#' @param beta Numeric vector of length J containing the infection-age spline
#'   coefficients.
#' @param subj_pre Named list of subject-level precomputations containing
#'   `V_i`, a numeric J x 366 cumulative-kernel matrix.
#' @param days Integer scalar in 1:365 giving the cumulative age endpoint.
#'   The default is 365.
#' @param control An `rsv_control` object controlling numerical bound checks
#'   and clipping behavior.
#'
#' @return Numeric scalar in [0, 1] giving the probability of first infection
#'   by age `days`.
#'
#' @details
#' For cumulative age endpoint \eqn{d}, the subject-specific cumulative hazard
#' is
#' \deqn{
#'   H_i(d) = \beta^\top v_i(d).
#' }
#' The probability of first infection by age \eqn{d} is therefore
#' \deqn{
#'   P(R_i \le d)
#'   =
#'   1 - \exp\{-H_i(d)\}.
#' }
#' The returned probability is clipped to [0, 1] according to the numerical
#' settings in `control`.
prob_infected_by_age <- function(beta,
                                 subj_pre,
                                 days = 365L,
                                 control = make_rsv_control()) {
  
  # Validate the age endpoint and cumulative-kernel dimensions.
  if (!is.numeric(days) || length(days) != 1L ||
      days < 1L || days > 365L)
    stop("days must be an integer in 1:365.")
  
  V_i <- subj_pre$V_i
  if (!is.matrix(V_i) || ncol(V_i) != 366L)
    stop("subj_pre$V_i must be J x 366 (from rsv_precompute_subject()).")
  
  if (length(beta) != nrow(V_i))
    stop("Length of beta must match nrow(V_i).")
  
  # Compute the cumulative hazard and first-infection probability.
  v_d <- V_i[, days + 1L]
  H_d <- sum(beta * v_d)
  
  Fbar_d <- exp(-H_d)
  
  prob <- 1 - Fbar_d
  
  # Enforce probability bounds using the configured numerical controls.
  .clip_to_range(
    prob,
    lower        = 0.0,
    upper        = 1.0,
    name         = "P(R_i <= days)",
    check_bounds = control$check_bounds,
    warn_on_clip = control$warn_on_clip
  )
}


#' Plot infection-age and healthcare-visit curve recovery
#'
#' Compares the estimated infection-age and healthcare-visit curves with the
#' known curves used to generate the simulated data. The plots are displayed
#' during interactive use and can optionally be saved as PNG files.
#'
#' @param curve_results Data frame containing the curve values to compare. It
#'   must include:
#' \describe{
#'   \item{\code{age_months}}{Numeric vector giving age in months.}
#'   \item{\code{w_true}}{Numeric vector containing the known simulation
#'     infection-age curve.}
#'   \item{\code{w_estimated}}{Numeric vector containing the estimated
#'     infection-age curve.}
#'   \item{\code{c_true}}{Numeric vector containing the known simulation
#'     healthcare-visit curve.}
#'   \item{\code{c_estimated}}{Numeric vector containing the estimated
#'     healthcare-visit curve.}
#' }
#' @param save_plots Logical scalar indicating whether to save the recovery
#'   plots as PNG files. The default is `TRUE`.
#' @param output_dir Character scalar giving the directory where saved figures
#'   are written. The default is `"outputs"`.
#'
#' @return `NULL`, returned invisibly.
#'
#' @details
#' During interactive use, the function displays separate comparisons for the
#' infection-age and healthcare-visit curves.
#'
#' When `save_plots = TRUE`, `output_dir` is created if necessary and the
#' figures are written as `infection_age_curve.png` and
#' `visit_age_curve.png`.
plot_curve_recovery <- function(
    curve_results,
    save_plots = TRUE,
    output_dir = "outputs"
) {
  required_columns <- c(
    "age_months",
    "w_true",
    "w_estimated",
    "c_true",
    "c_estimated"
  )
  
  if (
    !is.data.frame(curve_results) ||
    !all(required_columns %in% names(curve_results))
  ) {
    stop(
      "`curve_results` must contain the age, true-curve, and ",
      "estimated-curve columns."
    )
  }
  
  if (!is.logical(save_plots) || length(save_plots) != 1L || is.na(save_plots)) {
    stop("`save_plots` must be TRUE or FALSE.")
  }
  
  draw_w_plot <- function() {
    plot(
      x    = curve_results$age_months,
      y    = curve_results$w_estimated,
      type = "l",
      col  = "red",
      lwd  = 2,
      lty  = 2,
      ylim = range(
        curve_results$w_true,
        curve_results$w_estimated
      ),
      main = "Estimated vs. true infection-age curve, w(a)",
      xlab = "Age in months",
      ylab = "Age-specific infection weight"
    )
    
    lines(
      x   = curve_results$age_months,
      y   = curve_results$w_true,
      col = "black",
      lwd = 2,
      lty = 1
    )
    
    legend(
      "topright",
      legend = c("Estimated", "Ground truth"),
      col    = c("red", "black"),
      lwd    = 2,
      lty    = c(2, 1),
      bty    = "n"
    )
  }
  
  draw_c_plot <- function() {
    plot(
      x    = curve_results$age_months,
      y    = curve_results$c_estimated,
      type = "l",
      col  = "blue",
      lwd  = 2,
      lty  = 2,
      ylim = range(
        curve_results$c_true,
        curve_results$c_estimated
      ),
      main = "Estimated vs. true visit-age curve, c(a)",
      xlab = "Age in months",
      ylab = "Healthcare-visit weight"
    )
    
    lines(
      x   = curve_results$age_months,
      y   = curve_results$c_true,
      col = "black",
      lwd = 2,
      lty = 1
    )
    
    legend(
      "topright",
      legend = c("Estimated", "Ground truth"),
      col    = c("blue", "black"),
      lwd    = 2,
      lty    = c(2, 1),
      bty    = "n"
    )
  }
  
  # Display plots in interactive sessions.
  if (interactive()) {
    draw_w_plot()
    draw_c_plot()
  }
  
  # Save plots when requested, creating the output directory if needed.
  if (save_plots) {
    dir.create(
      output_dir,
      showWarnings = FALSE,
      recursive = TRUE
    )
    
    png(
      filename = file.path(output_dir, "infection_age_curve.png"),
      width    = 1080,
      height   = 672,
      res      = 120
    )
    draw_w_plot()
    dev.off()
    
    png(
      filename = file.path(output_dir, "visit_age_curve.png"),
      width    = 1080,
      height   = 672,
      res      = 120
    )
    draw_c_plot()
    dev.off()
    
    message(
      "Figures saved to\n",
      normalizePath(output_dir, mustWork = FALSE)
    )
  }
  
  invisible(NULL)
}
