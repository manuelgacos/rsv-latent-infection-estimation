#' Extract visit ages and summary columns from an rsv_data simulation
#'
#' @description
#' Converts an \code{rsv_data} object (a list of subject records) into a tidy,
#' one-row-per-subject data frame suitable for quick summaries and plots.
#'
#' The returned table always includes:
#' \itemize{
#'   \item subject identifiers (\code{i}, \code{id})
#'   \item birthday index (\code{birth_index}) plus calendar-style grouping columns
#'         (\code{birth_month}, \code{birth_decile})
#'   \item simulated visit outcome (\code{visit_age}) plus indicators
#'         (\code{visit_in_window}, \code{I_visit})
#'   \item visit-month-of-age (\code{visit_month}) for in-window visits
#' }
#'
#' Optionally (\code{compute_g=TRUE}), it also computes the model-implied visit
#' probability \eqn{g_i} using the same internal pipeline as
#' \code{run_empirical_checks()}:
#' \enumerate{
#'   \item build \code{Fbar_i} from \code{beta} and precomputed \code{V_i},
#'   \item build \code{pi} from \code{lambda_i} and the global \code{w} curve,
#'   \item build \code{Q} from \code{Fbar} and \code{pi},
#'   \item compute \code{U = S_day %*% Q} and then \code{g_i = eta^T U}.
#' }
#' When enabled, the table includes \code{g_i}, \code{I_minus_g},
#' \code{abs_I_minus_g}, and \code{g_decile}.
#'
#' @param data An \code{"rsv_data"} object with a \code{data$subjects} list.
#'   Each subject should include at least \code{birth_index} and \code{visit_age}.
#' @param model An \code{"rsv_model"} object. Needed for \code{days} and
#'   for computing \code{g_i} when \code{compute_g=TRUE}.
#' @param visit_age_max Integer upper bound for “in-window” visit ages.
#'   Defaults to \code{days}. Used to define \code{visit_in_window} and \code{I_visit}.
#' @param days Integer age window length. Defaults to \code{ncol(model$B_day)}.
#' @param compute_g Logical; if \code{TRUE}, compute per-subject \code{g_i}.
#'   Default \code{FALSE}.
#' @param beta Numeric vector of infection-age basis coefficients. Required if
#'   \code{compute_g=TRUE}.
#' @param eta Numeric vector of visit-age basis coefficients. Required if
#'   \code{compute_g=TRUE}.
#' @param control An \code{"rsv_control"} (or compatible list) passed to internal helpers
#'   used to compute \code{g_i}.
#' @param subj_pre_list Optional list of subject precomputations, one per subject,
#'   each containing \code{lambda_i} and \code{V_i}. REQUIRED if \code{compute_g=TRUE}.
#' @param month_len Numeric; length (in days) used to map day indices to a
#'   month-like bin via \code{ceiling(day / month_len)}. Default is \code{365/12}.
#' @param nbins Integer number of quantile bins for \code{g_decile} when
#'   \code{compute_g=TRUE}. Default is \code{10L}.
#'
#' @return A data.frame with one row per subject. Always includes:
#' \describe{
#'   \item{\code{i}}{Row index (1..n).}
#'   \item{\code{id}}{Subject id (if present in the record; otherwise NA).}
#'   \item{\code{birth_index}}{Subject birthday index used for calendar alignment.}
#'   \item{\code{birth_month}}{Month-like bin derived from \code{birth_index}.}
#'   \item{\code{birth_decile}}{Decile bin (1..10) of \code{birth_index}.}
#'   \item{\code{visit_age}}{Simulated visit day-of-age (NA if no visit).}
#'   \item{\code{visit_in_window}}{TRUE if \code{visit_age} is in \code{1..visit_age_max}.}
#'   \item{\code{I_visit}}{0/1 indicator equal to \code{as.integer(visit_in_window)}.}
#'   \item{\code{visit_month}}{Month-like bin of \code{visit_age} for in-window visits; NA otherwise.}
#' }
#' If \code{compute_g=TRUE}, also includes:
#' \describe{
#'   \item{\code{g_i}}{Model-implied probability of a visit within the evaluation window.}
#'   \item{\code{I_minus_g}}{Residual-like quantity \code{I_visit - g_i}.}
#'   \item{\code{abs_I_minus_g}}{Absolute deviation \code{abs(I_visit - g_i)}.}
#'   \item{\code{g_decile}}{Quantile bin (1..nbins) of \code{g_i}.}
#' }
#'
#' @examples
#' \dontrun{
#' # After simulating a dataset:
#' df <- extract_visit_ages(data_sim, model = sim$model, visit_age_max = 365)
#'
#' # With model-implied probabilities (requires subject precomputes):
#' subj_pre <- lapply(data_sim$subjects, rsv_precompute_subject, model = sim$model, control = make_rsv_control())
#' df2 <- extract_visit_ages(
#'   data_sim, model = sim$model, compute_g = TRUE,
#'   beta = beta_true, eta = eta_true, subj_pre_list = subj_pre
#' )
#' }
#'
#' @export
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
  
  # Helper: month index from day index (age day or birth day)
  day_to_month <- function(day, month_len, max_month = 12L) {
    ifelse(is.na(day), NA_integer_,
           pmin(max_month, pmax(1L, as.integer(ceiling(day / month_len)))))
  }
  
  # Helper: quantile bins (1..nbins). Robust to ties.
  quantile_bins <- function(x, nbins = 10L) {
    nbins <- as.integer(nbins)
    stopifnot(nbins >= 2L)
    
    out <- rep(NA_integer_, length(x))
    ok <- is.finite(x)
    if (!any(ok)) return(out)
    
    qs <- stats::quantile(x[ok], probs = seq(0, 1, length.out = nbins + 1L),
                          na.rm = TRUE, type = 7)
    
    # If breaks collapse (ties), fall back to rank-based binning
    if (any(diff(qs) <= 0)) {
      r <- rank(x[ok], ties.method = "average")
      out[ok] <- pmin(nbins, pmax(1L, ceiling(nbins * r / max(r))))
      return(out)
    }
    
    out[ok] <- as.integer(cut(x[ok], breaks = qs, include.lowest = TRUE, labels = FALSE))
    out
  }
  
  # ------------------------------------------------------------
  # Base extraction (always)
  # ------------------------------------------------------------
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
  
  # ------------------------------------------------------------
  # Optional: compute g_i (requires subj_pre_list + beta/eta)
  # ------------------------------------------------------------
  compute_g <- isTRUE(compute_g)
  if (compute_g) {
    if (is.null(subj_pre_list)) {
      stop("extract_visit_ages: compute_g=TRUE requires subj_pre_list (precomputed via rsv_precompute_subject).")
    }
    stopifnot(is.list(subj_pre_list), length(subj_pre_list) == n)
    stopifnot(is.numeric(beta), is.numeric(eta))
    
    # Global precompute for w_vec and c_vec
    glob <- rsv_precompute_global(beta = beta, eta = eta, model = model, control = control)
    w_vec <- glob$w_vec
    c_vec <- glob$c_vec
    
    # Defensive trim window supported by global objects
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


#' Probability of infection by a given age
#'
#' Computes the analytical probability
#'   P(R_i <= days) = 1 - exp{ - beta^T v_i(days) }
#' for a single subject under the RSV infection-age model.
#'
#' @param beta Numeric vector length J.
#' @param subj_pre Output of rsv_precompute_subject() for subject i.
#' @param days Integer in 1:365 (default 365).
#' @param control rsv_control object.
#'
#' @return Numeric scalar in [0,1]: probability of infection by age `days`.
#'
#' @export
prob_infected_by_age <- function(beta,
                                 subj_pre,
                                 days = 365L,
                                 control = make_rsv_control()) {
  
  # ---- checks ----
  if (!is.numeric(days) || length(days) != 1L ||
      days < 1L || days > 365L)
    stop("days must be an integer in 1:365.")
  
  V_i <- subj_pre$V_i
  if (!is.matrix(V_i) || ncol(V_i) != 366L)
    stop("subj_pre$V_i must be J x 366 (from rsv_precompute_subject()).")
  
  if (length(beta) != nrow(V_i))
    stop("Length of beta must match nrow(V_i).")
  
  # ---- integrated hazard at 'days' ----
  v_d <- V_i[, days + 1L]      # v_i(days)
  H_d <- sum(beta * v_d)
  
  # ---- survival and probability ----
  Fbar_d <- exp(-H_d)
  
  prob <- 1 - Fbar_d
  
  # numerical guard
  .clip_to_range(
    prob,
    lower        = 0.0,
    upper        = 1.0,
    name         = "P(R_i <= days)",
    check_bounds = control$check_bounds,
    warn_on_clip = control$warn_on_clip
  )
}


# Scale factor s so that P_i(infected by D) = p_target for one subject i
scale_w_for_infection_prob_subject <- function(p_target,
                                               lambda_i,
                                               w0,
                                               D = 365L) {
  stopifnot(is.numeric(p_target), length(p_target) == 1L,
            is.finite(p_target), p_target > 0, p_target < 1)
  stopifnot(is.numeric(lambda_i), length(lambda_i) >= D)
  stopifnot(is.numeric(w0),       length(w0)       >= D)
  
  A <- sum(lambda_i[1:D] * w0[1:D])
  if (!is.finite(A) || A <= 0) {
    stop("A = sum(lambda_i * w0) must be positive and finite; check lambda_i/w0.")
  }
  
  s <- -log1p(-p_target) / A  # stable version of -log(1 - p)
  s
}


# Optimized, fully vectorized version
scale_c_for_visit_prob_subject <- function(p_target,
                                           lambda_i,
                                           w,
                                           c0,
                                           D = 365L) {
  stopifnot(length(lambda_i) >= D,
            length(w) >= D,
            length(c0) >= D,
            is.numeric(p_target), length(p_target) == 1L,
            is.finite(p_target), p_target > 0, p_target < 1)
  
  hazard <- lambda_i[1:D] * w[1:D]
  cumhazard <- c(0, cumsum(hazard))[1:D]  # lagged version
  Fbar <- exp(-cumhazard)
  
  pi <- hazard * Fbar
  B <- sum(pi * c0[1:D])
  
  if (B <= 0 || !is.finite(B))
    stop("B must be positive and finite. Check inputs.")
  
  s <- p_target / B
  s
}


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
      ylab = "Probability of visit given infection"
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
  
  # Display the plots during interactive use.
  if (interactive()) {
    draw_w_plot()
    draw_c_plot()
  }
  
  # Save the plots for reproducible terminal and interactive runs.
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