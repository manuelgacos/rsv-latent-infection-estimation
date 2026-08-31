# Probability core
#
# Implements the subject-level precomputations and probability calculations
# connecting seasonal RSV circulation, first-infection timing, and healthcare visits.

#' Compute subject-level RSV precomputations
#'
#' Computes the subject-specific quantities that depend on birth timing and the
#' model but not on the infection-age or healthcare-visit spline coefficients.
#' These quantities are reused across likelihood and estimation evaluations.
#'
#' @param subject_i A subject record containing `birth_index`, the 1-based
#'   calendar index used to align the seasonal RSV circulation curve with the
#'   subject's age.
#' @param model An `rsv_model` object containing the daily infection-age basis
#'   `B_day`, integrated basis `phi`, and seasonal RSV circulation curve
#'   `lambda_global`.
#' @param control An `rsv_control` object created by `make_rsv_control()`,
#'   controlling calendar-index bound checks and clipping behavior.
#'
#' @return A named list of subject-level precomputations containing:
#' \describe{
#'   \item{\code{lambda_i}}{Numeric vector of length 365 containing the
#'     subject-specific RSV circulation curve over ages 1 through 365 days.}
#'   \item{\code{V_i}}{Numeric \code{J x 366} matrix containing the cumulative
#'     subject-level kernels. Column 1 represents day zero, and column
#'     \code{d + 1} represents \eqn{v_i(d)} for days 1 through 365.}
#' }
#'
#' @details
#' The subject-specific circulation curve is obtained by shifting
#' \eqn{\lambda(t)} according to `subject_i$birth_index`. The cumulative kernel
#' matrix is then computed from `lambda_i` and `model$phi`. Neither returned
#' quantity depends on \eqn{\beta} or \eqn{\eta}.
rsv_precompute_subject <- function(subject_i,
                                   model,
                                   control = make_rsv_control()) {
  if (!is.list(subject_i)) {
    stop("`subject_i` must be a list (typically from make_subject()).")
  }
  if (!inherits(model, "rsv_model")) {
    stop("`model` must be an object of class 'rsv_model' (from make_rsv_model()).")
  }
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("`control` must be an 'rsv_control' object (from make_rsv_control()).")
  }
  if (is.null(model$phi)) {
    stop("model$phi is missing; it must be provided in the rsv_model object.")
  }
  if (is.null(model$lambda_global)) {
    stop("model$lambda_global is missing; it must be provided in the rsv_model object.")
  }
  
  phi <- model$phi
  
  if (!is.matrix(phi) || !is.numeric(phi)) {
    stop("model$phi must be a numeric matrix.")
  }
  if (ncol(phi) != 365L) {
    stop(sprintf("model$phi must have 365 columns (found %d).", ncol(phi)))
  }
  
  # Construct the subject-specific circulation curve.
  lambda_i <- lambda_shift_i(subject_i, model, control)
  
  if (!is.numeric(lambda_i)) {
    stop("lambda_shift_i() must return a numeric vector.")
  }
  if (length(lambda_i) != ncol(phi)) {
    stop(sprintf(
      "Length of lambda_i (%d) must match ncol(model$phi) (%d).",
      length(lambda_i), ncol(phi)
    ))
  }
  
  # Construct the subject-specific cumulative kernel.
  V_i <- v_i(lambda_i, phi)
  
  if (!is.matrix(V_i) || !is.numeric(V_i)) {
    stop("v_i() must return a numeric matrix.")
  }
  if (ncol(V_i) != 366L) {
    stop(sprintf(
      "V_i must have 366 columns (v_i(0),...,v_i(365)); found %d.",
      ncol(V_i)
    ))
  }
  
  list(
    lambda_i = lambda_i,
    V_i      = V_i
  )
}


#' Construct the subject-specific RSV circulation curve
#'
#' Aligns the seasonal RSV circulation curve \eqn{\lambda(t)} with a subject's
#' age using the subject's birth index.
#'
#' @param subject_i A subject record containing `birth_index`, the 1-based
#'   calendar index corresponding to the subject's birth date.
#' @param model An `rsv_model` object containing the seasonal RSV circulation
#'   curve `lambda_global` and the daily age basis `B_day`.
#' @param control An `rsv_control` object controlling calendar-index bound
#'   checks and clipping behavior.
#'
#' @return Numeric vector of length 365 containing the subject-specific RSV
#'   circulation values over ages 1 through 365 days.
#'
#' @details
#' For birth index \eqn{B_i} and age \eqn{a}, the subject-specific circulation
#' curve is
#' \deqn{
#'   \lambda_i(a) = \lambda(B_i + a).
#' }
#' Thus, ages 1 through 365 correspond to calendar indices
#' `B_i + 1` through `B_i + 365`.
#'
#' If a shifted calendar index falls outside `lambda_global`, strict mode
#' produces an error. Otherwise, the index is clipped to the available
#' calendar range, with an optional warning.
lambda_shift_i <- function(subject_i, model, control) {
  
  B_i <- subject_i$birth_index
  lambda_global <- model$lambda_global
  
  # Infer the modeled age window from the daily basis.
  age_len  <- ncol(model$B_day)
  age_grid <- seq_len(age_len)
  
  # Map ages 1:age_len to calendar indices B_i + 1 through B_i + age_len.
  cal_idx <- B_i + age_grid
  
  if (control$check_bounds) {
    L <- length(lambda_global)
    if (any(cal_idx < 1L | cal_idx > L)) {
      stop(sprintf(
        "lambda_shift_i: calendar indices out of bounds for subject with birth_index = %d.",
        B_i
      ))
    }
  } else {
    # In lenient mode, clip out-of-range calendar indices.
    L <- length(lambda_global)
    too_low  <- cal_idx < 1L
    too_high <- cal_idx > L
    if (any(too_low | too_high)) {
      if (control$warn_on_clip) {
        warning("lambda_shift_i: calendar indices out of bounds; clipping applied.")
      }
      cal_idx <- pmin(pmax(cal_idx, 1L), L)
    }
  }
  
  lambda_global[cal_idx]
}


#' Compute the subject-specific cumulative kernel
#'
#' Computes the cumulative kernel \eqn{v_i(d)} from the subject-specific RSV
#' circulation curve and the integrated infection-age spline basis.
#'
#' @param lambda_shift_i Numeric vector of length 365 containing the
#'   subject-specific RSV circulation values over ages 1 through 365 days.
#' @param phi Numeric `J x 365` matrix containing the integrated infection-age
#'   spline basis, where column \eqn{m} corresponds to age day \eqn{m}.
#'
#' @return Numeric `J x 366` matrix containing the cumulative kernels.
#'   Column 1 represents \eqn{v_i(0) = 0}, and column `d + 1` represents
#'   \eqn{v_i(d)} for days 1 through 365.
#'
#' @details
#' For age day \eqn{d},
#' \deqn{
#'   v_i(d)
#'   =
#'   \sum_{m=1}^{d} \lambda_i(m)\phi(m).
#' }
#' The additional first column preserves the day-zero value needed by
#' downstream survival and likelihood calculations.
v_i <- function(lambda_shift_i, phi) {
  if (!is.numeric(lambda_shift_i) || length(lambda_shift_i) != 365L) {
    stop("lambda_shift_i must be a numeric vector of length 365.")
  }
  if (!is.matrix(phi) || ncol(phi) != 365L) {
    stop("phi must be a numeric J x 365 matrix (columns are phi[, m]).")
  }
  J <- nrow(phi)
  
  # Scale each integrated basis column by the subject-specific circulation.
  A <- sweep(phi, 2L, lambda_shift_i, "*")
  
  # Accumulate v_i(d), with column 1 reserved for v_i(0) = 0.
  V <- matrix(0.0, nrow = J, ncol = 366L)
  for (m in 1:365) {
    V[, m + 1L] <- V[, m] + A[, m]
  }
  
  rownames(V) <- rownames(phi)
  colnames(V) <- c("d000", sprintf("d%03d", 1:365))
  V
}


#' Compute global RSV precomputations
#'
#' Computes the daily infection-age and healthcare-visit curves for the current
#' spline coefficients. These parameter-dependent quantities are reused across
#' subjects during likelihood and estimation calculations.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param model An `rsv_model` object containing the daily spline bases
#'   `B_day` and `S_day`.
#' @param control An `rsv_control` object controlling nonnegativity checks and
#'   numerical stabilization.
#'
#' @return A named list of global precomputations containing:
#' \describe{
#'   \item{\code{w_vec}}{Numeric vector of length 365 containing the
#'     infection-age curve \eqn{w(a)} over ages 1 through 365 days.}
#'   \item{\code{c_vec}}{Numeric vector of length 365 containing the
#'     healthcare-visit curve \eqn{c(a)} over ages 1 through 365 days.}
#' }
#'
#' @details
#' The curves are obtained by projecting \eqn{\beta} and \eqn{\eta} onto their
#' corresponding daily spline bases. Numerical safeguards are applied through
#' `control` before the curves are returned.
rsv_precompute_global <- function(beta,
                                  eta,
                                  model,
                                  control = make_rsv_control()) {
  if (!inherits(model, "rsv_model")) {
    stop("`model` must be an object of class 'rsv_model' (from make_rsv_model()).")
  }
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("`control` must be an 'rsv_control' object (from make_rsv_control()).")
  }
  
  B_day <- model$B_day
  S_day <- model$S_day
  
  if (!is.numeric(B_day) || !(is.matrix(B_day) || is.array(B_day))) {
    stop("model$B_day must be a numeric matrix or array.")
  }
  if (!is.numeric(S_day) || !(is.matrix(S_day) || is.array(S_day))) {
    stop("model$S_day must be a numeric matrix or array.")
  }
  if (ncol(B_day) != 365L) {
    stop(sprintf("model$B_day must have 365 columns (found %d).", ncol(B_day)))
  }
  if (ncol(S_day) != 365L) {
    stop(sprintf("model$S_day must have 365 columns (found %d).", ncol(S_day)))
  }
  if (length(beta) != nrow(B_day)) {
    stop(sprintf(
      "Length of beta (%d) must match nrow(model$B_day) (%d).",
      length(beta), nrow(B_day)
    ))
  }
  if (length(eta) != nrow(S_day)) {
    stop(sprintf(
      "Length of eta (%d) must match nrow(model$S_day) (%d).",
      length(eta), nrow(S_day)
    ))
  }
  
  # Evaluate the daily infection-age and healthcare-visit curves.
  w_vec_raw <- w_day(beta, B_day)
  c_vec_raw <- c_day(eta,  S_day)
  
  # Enforce the model's nonnegativity constraints.
  w_vec <- .enforce_nonneg(
    x       = w_vec_raw,
    name    = "w_vec",
    control = control
  )
  
  c_vec <- .enforce_nonneg(
    x       = c_vec_raw,
    name    = "c_vec",
    control = control
  )
  
  # Apply a positive floor before downstream logarithms.
  w_vec <- .stabilize_for_log(
    x            = w_vec,
    name         = "w_vec",
    eps_log      = control$eps_log,
    check_bounds = control$check_bounds,
    warn_on_clip = control$warn_on_clip
  )
  
  c_vec <- .stabilize_for_log(
    x            = c_vec,
    name         = "c_vec",
    eps_log      = control$eps_log,
    check_bounds = control$check_bounds,
    warn_on_clip = control$warn_on_clip
  )
  
  list(
    w_vec = w_vec,
    c_vec = c_vec
  )
}


#' Evaluate the infection-age curve
#'
#' Evaluates the infection-age curve \eqn{w(a)} on the daily age grid from the
#' infection-age spline coefficients and basis matrix.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param B_day Numeric `J x 365` matrix containing the daily infection-age
#'   spline basis, where column \eqn{a} is the basis vector \eqn{b(a)}.
#'
#' @return Numeric vector of length 365 containing the infection-age curve
#'   evaluated over ages 1 through 365 days.
#'
#' @details
#' For age day \eqn{a},
#' \deqn{
#'   w(a) = \beta^\top b(a).
#' }
w_day <- function(beta, B_day) {
  if (!is.numeric(beta) || !is.numeric(B_day)) {
    stop("beta and B_day must be numeric.")
  }
  if (length(beta) != nrow(B_day)) {
    stop("Length of beta must match number of rows in B_day.")
  }
  # Evaluate beta^T b(a) for each age day.
  as.numeric(crossprod(B_day, beta))
}


#' Evaluate the healthcare-visit curve
#'
#' Evaluates the healthcare-visit curve \eqn{c(a)} on the daily age grid from
#' the healthcare-visit spline coefficients and basis matrix.
#'
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param S_day Numeric `K x 365` matrix containing the daily healthcare-visit
#'   spline basis, where column \eqn{a} is the basis vector \eqn{s(a)}.
#'
#' @return Numeric vector of length 365 containing the healthcare-visit curve
#'   evaluated over ages 1 through 365 days.
#'
#' @details
#' For age day \eqn{a},
#' \deqn{
#'   c(a) = \eta^\top s(a).
#' }
c_day <- function(eta, S_day) {
  if (!is.numeric(eta) || !is.numeric(S_day)) {
    stop("eta and S_day must be numeric.")
  }
  if (length(eta) != nrow(S_day)) {
    stop("Length of eta must match number of rows in S_day.")
  }
  # Evaluate eta^T s(a) for each age day.
  as.numeric(crossprod(S_day, eta))
}


#' Enforce nonnegativity for numeric values
#'
#' Replaces negative entries with zero while optionally distinguishing small
#' numerical deviations from larger violations using `control$tol_clip`.
#'
#' @param x Numeric vector expected to be nonnegative.
#' @param name Character scalar used in warning and error messages.
#' @param control An `rsv_control` object controlling the clipping tolerance,
#'   strict bound checks, and warning behavior.
#'
#' @return Numeric vector of the same length as `x`, with negative entries
#'   replaced by zero.
#'
#' @details
#' Values below `-control$tol_clip` produce an error when strict bound checking
#' is enabled. Otherwise, all negative values are clipped to zero; larger
#' violations always produce a warning, while small numerical deviations warn
#' only when `control$warn_on_clip` is enabled.
.enforce_nonneg <- function(x,
                            name    = "value",
                            control = make_rsv_control()) {
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric in .enforce_nonneg().", name))
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` contains non-finite values (NA, NaN, or Inf).", name))
  }
  
  # Use defaults when individual control fields are unavailable.
  tol_clip     <- if (!is.null(control$tol_clip))     control$tol_clip     else 1e-12
  check_bounds <- if (!is.null(control$check_bounds)) control$check_bounds else FALSE
  warn_on_clip <- if (!is.null(control$warn_on_clip)) control$warn_on_clip else TRUE
  
  if (!is.numeric(tol_clip) || length(tol_clip) != 1L || tol_clip < 0 || !is.finite(tol_clip)) {
    stop("control$tol_clip must be a single nonnegative finite numeric value.")
  }
  
  neg_idx <- which(x < 0)
  if (length(neg_idx) == 0L) {
    return(x)
  }
  
  x_neg <- x[neg_idx]
  
  # Separate numerical drift from larger nonnegativity violations.
  tiny_idx    <- neg_idx[x_neg >= -tol_clip]
  serious_idx <- neg_idx[x_neg < -tol_clip]
  
  # In strict mode, larger violations are treated as errors.
  if (length(serious_idx) > 0L && isTRUE(check_bounds)) {
    min_val <- min(x[serious_idx])
    stop(sprintf(
      "%s has %d values < -tol_clip (min = %g) violating nonnegativity.",
      name, length(serious_idx), min_val
    ))
  }
  
  if (length(tiny_idx) > 0L) {
    x[tiny_idx] <- 0
    if (isTRUE(warn_on_clip)) {
      warning(sprintf(
        "%s had %d small negative values in [-tol_clip, 0); clipped to 0.",
        name, length(tiny_idx)
      ))
    }
  }
  
  if (length(serious_idx) > 0L) {
    x[serious_idx] <- 0
    # Larger violations always warn when strict checking is disabled.
    min_val <- min(x_neg[x_neg < -tol_clip])
    warning(sprintf(
      "%s had %d values < -tol_clip (min = %g); clipped to 0. This may indicate an issue with the model or parameters.",
      name, length(serious_idx), min_val
    ))
  }
  
  x
}


#' Stabilize nonnegative values for logarithms
#'
#' Applies a strictly positive lower bound to numeric values before they are
#' used inside logarithms.
#'
#' @param x Numeric vector expected to be nonnegative.
#' @param name Character scalar used in warning and error messages.
#' @param eps_log Positive numeric scalar giving the lower bound applied to
#'   values below `eps_log`.
#' @param check_bounds Logical scalar controlling whether negative values
#'   produce an error.
#' @param warn_on_clip Logical scalar controlling warnings when values are
#'   raised to `eps_log`.
#'
#' @return Numeric vector of the same length as `x`, with all entries at least
#'   `eps_log`.
#'
#' @details
#' Negative values produce an error when `check_bounds = TRUE`; otherwise they
#' are raised to `eps_log`. Nonnegative values below `eps_log` are also raised
#' to the lower bound.
.stabilize_for_log <- function(x,
                               name,
                               eps_log,
                               check_bounds = FALSE,
                               warn_on_clip = TRUE) {
  
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric in .stabilize_for_log().", name))
  }
  if (!is.numeric(eps_log) || length(eps_log) != 1L ||
      !is.finite(eps_log) || eps_log <= 0) {
    stop("`eps_log` must be a single positive finite numeric value in .stabilize_for_log().")
  }
  
  neg_idx <- which(x < 0)
  if (length(neg_idx) > 0L) {
    if (check_bounds) {
      stop(sprintf(
        "%s has %d negative values in .stabilize_for_log(); input is expected to be nonnegative.",
        name, length(neg_idx)
      ))
    } else {
      # In lenient mode, repair negative values at the logarithmic floor.
      if (warn_on_clip) {
        warning(sprintf(
          "%s has %d negative values; stabilizing them to eps_log = %g in .stabilize_for_log().",
          name, length(neg_idx), eps_log
        ))
      }
      x[neg_idx] <- eps_log
    }
  }
  
  # Raise small nonnegative values to the strictly positive logarithmic floor.
  small_idx <- which(x >= 0 & x < eps_log)
  if (length(small_idx) > 0L) {
    if (warn_on_clip) {
      warning(sprintf(
        "%s has %d values in [0, eps_log); stabilizing them to eps_log = %g for log-safety.",
        name, length(small_idx), eps_log
      ))
    }
    x[small_idx] <- eps_log
  }
  
  x
}


#' Compute the subject-specific cumulative hazard
#'
#' Computes the cumulative hazard \eqn{H_i(d)} from the infection-age spline
#' coefficients and the subject-specific cumulative kernel.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param V_i Numeric `J x 366` matrix containing the cumulative subject-level
#'   kernels. Column 1 represents \eqn{v_i(0) = 0}, and column `d + 1`
#'   represents \eqn{v_i(d)} for days 1 through 365.
#' @param include_day0 Logical scalar controlling whether the returned vector
#'   includes \eqn{H_i(0) = 0}.
#'
#' @return Numeric vector of length 365 containing \eqn{H_i(d)} for days
#'   1 through 365, or length 366 with day zero prepended when
#'   `include_day0 = TRUE`.
#'
#' @details
#' For age day \eqn{d},
#' \deqn{
#'   H_i(d) = \beta^\top v_i(d).
#' }
H_i <- function(beta, V_i, include_day0 = FALSE) {
  if (!is.numeric(beta) || !is.numeric(V_i)) {
    stop("beta and V_i must be numeric.")
  }
  if (length(beta) != nrow(V_i)) {
    stop("Length of beta must match number of rows in V_i.")
  }
  if (ncol(V_i) != 366) {
    stop("V_i must have 366 columns: v_i(0), v_i(1), ..., v_i(365).")
  }
  # Skip the day-zero column and evaluate beta^T v_i(d) for days 1:365.
  H <- as.numeric(crossprod(V_i[, 2:366, drop = FALSE], beta))
  if (isTRUE(include_day0)) {
    H <- c(0, H)
  }
  H
}


#' Convert cumulative hazards to survival values
#'
#' Computes subject-specific survival values from cumulative hazards and
#' applies numerical clipping to preserve valid survival bounds.
#'
#' @param H Numeric vector of length 365 containing \eqn{H_i(d)} for days
#'   1 through 365, or length 366 including \eqn{H_i(0) = 0} when
#'   `include_day0 = TRUE`.
#' @param include_day0 Logical scalar indicating whether `H` includes day zero.
#' @param check_bounds Logical scalar controlling strict bound checks during
#'   numerical clipping.
#' @param warn_on_clip Logical scalar controlling warnings when clipping occurs.
#'
#' @return Numeric vector of the same length as `H` containing the
#'   subject-specific survival values.
#'
#' @details
#' Survival is computed as
#' \deqn{
#'   \bar F_i(d) = \exp\{-H_i(d)\}.
#' }
#' Values are clipped to \eqn{(0, 1]} using machine precision as the positive
#' lower bound.
Fbar_from_H <- function(H,
                        include_day0 = FALSE,
                        check_bounds = FALSE,
                        warn_on_clip = TRUE) {
  
  if (!is.numeric(H))
    stop("H must be numeric.")
  
  len_H <- length(H)
  
  if (include_day0) {
    # Day-zero input must begin with H_i(0) = 0.
    if (len_H != 366)
      stop("If include_day0=TRUE, H must have length 366 (H_i(0)...H_i(365)).")
    
    Fbar <- exp(-H)
    
  } else {
    if (len_H != 365)
      stop("If include_day0=FALSE, H must have length 365 (H_i(1)...H_i(365)).")
    
    Fbar <- exp(-H)
  }
  
  # Keep survival values in (0, 1] and avoid numerical underflow to zero.
  eps <- .Machine$double.eps
  
  Fbar <- .clip_to_range(
    Fbar,
    lower        = eps,
    upper        = 1.0,
    name         = "Fbar_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  Fbar
}


#' Compute the subject-specific survival curve
#'
#' Computes the probability that a subject remains uninfected through each age
#' day from the infection-age spline coefficients and subject-specific
#' cumulative kernel.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param V_i Numeric `J x 366` matrix containing the cumulative subject-level
#'   kernels. Column 1 represents \eqn{v_i(0) = 0}, and column `d + 1`
#'   represents \eqn{v_i(d)} for days 1 through 365.
#' @param include_day0 Logical scalar controlling whether the returned vector
#'   includes \eqn{\bar F_i(0) = 1}.
#' @param check_bounds Logical scalar controlling strict bound checks during
#'   survival-value clipping.
#' @param warn_on_clip Logical scalar controlling warnings when clipping is
#'   applied.
#'
#' @return Numeric vector of length 365 containing \eqn{\bar F_i(d)} for days
#'   1 through 365, or length 366 with day zero prepended when
#'   `include_day0 = TRUE`.
#'
#' @details
#' The survival curve is
#' \deqn{
#'   \bar F_i(d)
#'   =
#'   \exp\{-H_i(d)\}
#'   =
#'   \exp\{-\beta^\top v_i(d)\}.
#' }
#' The cumulative hazards are computed by `H_i()`, and numerical clipping is
#' applied by `Fbar_from_H()`.
Fbar_i <- function(beta,
                   V_i,
                   include_day0 = FALSE,
                   check_bounds = FALSE,
                   warn_on_clip = TRUE) {
  if (!is.numeric(beta) || !is.numeric(V_i)) {
    stop("beta and V_i must be numeric.")
  }
  if (length(beta) != nrow(V_i)) {
    stop("Length of beta must match number of rows in V_i.")
  }
  
  H <- H_i(beta, V_i, include_day0 = include_day0)
  
  Fbar_from_H(
    H,
    include_day0 = include_day0,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
}


#' Compute daily first-infection probabilities
#'
#' Computes the subject-specific probability of first infection on each age day,
#' conditional on remaining uninfected through the previous day.
#'
#' @param lambda_shift_i Numeric vector of length 365 containing the
#'   subject-specific RSV circulation values over ages 1 through 365 days.
#' @param w Numeric vector of length 365 containing the infection-age curve
#'   \eqn{w(a)} over ages 1 through 365 days.
#' @param check_bounds Logical scalar controlling strict probability-bound
#'   checks during numerical clipping.
#' @param warn_on_clip Logical scalar controlling warnings when clipping is
#'   applied.
#'
#' @return Numeric vector of length 365 containing the conditional
#'   first-infection probabilities \eqn{\pi_i(a)}.
#'
#' @details
#' For age day \eqn{a},
#' \deqn{
#'   \pi_i(a)
#'   =
#'   1 - \exp\{-\lambda_i(a)w(a)\}.
#' }
#' The probabilities are clipped to \eqn{[0, 1]} to guard against numerical
#' drift outside the valid probability range.
pi_from_w <- function(lambda_shift_i,
                      w,
                      check_bounds = FALSE,
                      warn_on_clip = TRUE) {
  if (!is.numeric(lambda_shift_i) || length(lambda_shift_i) != 365L)
    stop("lambda_shift_i must be numeric length 365.")
  if (!is.numeric(w) || length(w) != 365L)
    stop("w must be numeric length 365.")
  
  # Use expm1() for stable evaluation of 1 - exp(-x).
  x  <- lambda_shift_i * w
  pi <- -expm1(-x)
  
  # Guard against numerical drift outside the probability range.
  pi <- .clip_to_range(
    pi,
    lower        = 0.0,
    upper        = 1.0,
    name         = "pi_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  pi
}


#' Compute the subject-specific first-infection mass
#'
#' Computes the unconditional probability mass of first infection on each age
#' day from the subject-specific survival curve and conditional first-infection
#' probabilities.
#'
#' @param Fbar_i Numeric vector of length 366 containing the subject-specific
#'   survival curve from day zero through day 365, where element `d + 1`
#'   represents \eqn{\bar F_i(d)}.
#' @param pi_i Numeric vector of length 365 containing the conditional
#'   first-infection probabilities \eqn{\pi_i(a)} for ages 1 through 365 days.
#' @param check_bounds Logical scalar controlling strict probability-bound
#'   checks during numerical clipping.
#' @param warn_on_clip Logical scalar controlling warnings when clipping is
#'   applied.
#'
#' @return Numeric vector of length 365 containing the first-infection
#'   probability masses \eqn{Q_i(a)}.
#'
#' @details
#' For age day \eqn{a},
#' \deqn{
#'   Q_i(a)
#'   =
#'   \bar F_i(a - 1)\pi_i(a).
#' }
#' The resulting values are clipped to \eqn{[0, 1]} to guard against numerical
#' drift outside the valid probability range.
Q_i <- function(Fbar_i,
                pi_i,
                check_bounds = FALSE,
                warn_on_clip = TRUE) {
  if (!is.numeric(Fbar_i) || length(Fbar_i) != 366L) {
    stop("Fbar_i must be a numeric vector of length 366 (d = 0..365).")
  }
  if (!is.numeric(pi_i) || length(pi_i) != 365L) {
    stop("pi_i must be a numeric vector of length 365 (m = 1..365).")
  }
  
  # Enforce the theoretical bounds for survival and conditional probabilities.
  eps <- .Machine$double.eps
  
  Fbar_i <- .clip_to_range(
    Fbar_i,
    lower        = eps,
    upper        = 1.0,
    name         = "Fbar_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  pi_i <- .clip_to_range(
    pi_i,
    lower        = 0.0,
    upper        = 1.0,
    name         = "pi_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  # Element m of Fbar_i represents survival through day m - 1.
  Q <- Fbar_i[1:365] * pi_i
  
  # Guard against numerical drift outside the probability range.
  Q <- .clip_to_range(
    Q,
    lower        = 0.0,
    upper        = 1.0,
    name         = "Q_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  Q
}


#' Compute the subject-specific healthcare-visit kernel
#'
#' Projects the subject-specific first-infection probability mass onto the
#' healthcare-visit spline basis.
#'
#' @param Q_i Numeric vector of length 365 containing the first-infection
#'   probability masses \eqn{Q_i(a)} over ages 1 through 365 days.
#' @param S_day Numeric `K x 365` matrix containing the daily healthcare-visit
#'   spline basis, where column \eqn{a} is the basis vector \eqn{s(a)}.
#' @param check_bounds Logical scalar controlling strict bound checks when
#'   validating `Q_i`.
#' @param warn_on_clip Logical scalar controlling warnings when clipping is
#'   applied to `Q_i`.
#'
#' @return Numeric vector of length `K` containing the subject-specific
#'   healthcare-visit kernel \eqn{U_i}.
#'
#' @details
#' The kernel is
#' \deqn{
#'   U_i
#'   =
#'   \sum_{a=1}^{365} s(a)Q_i(a).
#' }
#' Equivalently, \eqn{U_i = S Q_i}. Small negative values introduced by
#' numerical error are set to zero.
U_i <- function(Q_i,
                S_day,
                check_bounds = FALSE,
                warn_on_clip = TRUE) {
  if (!is.numeric(Q_i) || length(Q_i) != 365L) {
    stop("Q_i must be a numeric vector of length 365.")
  }
  
  # Accept either a base R matrix or a Matrix-class object.
  isValidMatrix <-
    (is.matrix(S_day) ||
       inherits(S_day, "Matrix")) &&
    is.numeric(S_day)
  
  if (!isValidMatrix) {
    stop("S_day must be a numeric matrix (base R or Matrix package).")
  }
  
  if (ncol(S_day) != 365L) {
    stop("S_day must have 365 columns (one for each s(m)).")
  }
  
  # Enforce the probability bounds for Q_i before projection.
  Q_i <- .clip_to_range(
    Q_i,
    lower        = 0.0,
    upper        = 1.0,
    name         = "Q_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  U <- as.vector(S_day %*% Q_i)
  
  # Remove tiny negative values introduced by numerical error.
  U[U < 0 & U > -1e-12] <- 0
  
  U
}


#' Compute the subject-specific no-visit probability
#'
#' Computes the probability that a subject has no healthcare visit during the
#' first year of life from the subject-specific healthcare-visit kernel and
#' healthcare-visit spline coefficients.
#'
#' @param U_i Numeric vector of length `K` containing the subject-specific
#'   healthcare-visit kernel.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param eps_tau Nonnegative numeric scalar giving the lower bound applied to
#'   the no-visit probability. Set to zero to disable the positive floor.
#' @param check_bounds Logical scalar controlling strict bound checks during
#'   numerical clipping.
#' @param warn_on_clip Logical scalar controlling warnings when clipping is
#'   applied.
#'
#' @return Numeric scalar containing the no-visit probability \eqn{\tau_i},
#'   bounded to the interval \eqn{[\mathrm{eps\_tau}, 1]}.
#'
#' @details
#' The no-visit probability is
#' \deqn{
#'   \tau_i(\beta, \eta)
#'   =
#'   1 - \eta^\top U_i(\beta).
#' }
#' A lower bound of `eps_tau` is applied to avoid numerical instability when
#' the probability is close to zero.
tau_i <- function(U_i,
                  eta,
                  eps_tau = 1e-12,
                  check_bounds = FALSE,
                  warn_on_clip = TRUE) {
  if (!is.numeric(U_i) || !is.numeric(eta))
    stop("U_i and eta must be numeric.")
  if (length(U_i) != length(eta))
    stop("U_i and eta must have the same length (K).")
  if (!is.finite(eps_tau) || eps_tau < 0)
    stop("eps_tau must be a nonnegative finite number.")
  
  dot <- sum(eta * U_i)
  if (!is.finite(dot))
    stop("Non-finite dot product: check inputs.")
  
  tau <- 1 - dot
  
  # Use eps_tau rather than machine epsilon to avoid unstable 
  # near-zero no-visit probabilities.
  tau <- .clip_to_range(
    tau,
    lower        = eps_tau,
    upper        = 1.0,
    name         = "tau_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  tau
}


#' Clip numeric values to a bounded interval
#'
#' Clips numeric values to a specified interval, with optional checks for
#' violations beyond a numerical tolerance.
#'
#' @param x Numeric vector containing the values to clip.
#' @param lower,upper Numeric scalars defining the clipping interval, with
#'   `lower <= upper`.
#' @param name Character scalar used in warning and error messages.
#' @param check_bounds Logical scalar controlling whether violations beyond
#'   `tol` produce an error.
#' @param warn_on_clip Logical scalar controlling warnings when clipping is
#'   applied.
#' @param tol Nonnegative numeric scalar distinguishing small numerical
#'   deviations from larger bound violations.
#'
#' @return Numeric vector of the same length as `x`, with values clipped to
#'   the interval [`lower`, `upper`].
#'
#' @details
#' Values outside the interval by more than `tol` produce an error when
#' `check_bounds = TRUE`. Otherwise, out-of-range values are clipped to the
#' specified bounds.
.clip_to_range <- function(x, lower, upper,
                           name = "value",
                           check_bounds = FALSE,
                           warn_on_clip = TRUE,
                           tol = 1e-12) {
  
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric in .clip_to_range().", name))
  }
  if (!is.numeric(lower) || !is.numeric(upper) ||
      length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper) ||
      lower > upper) {
    stop("Invalid bounds: 'lower' and 'upper' must be finite scalars with lower <= upper.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || tol < 0 || !is.finite(tol)) {
    stop("Argument 'tol' must be a nonnegative finite numeric scalar.")
  }
  
  below <- x < lower
  above <- x > upper
  any_oob <- any(below | above)
  
  if (!any_oob) {
    return(x)
  }
  
  # Distinguish numerical drift from violations beyond tolerance.
  far_below <- x < (lower - tol)
  far_above <- x > (upper + tol)
  n_far <- sum(far_below | far_above)
  
  if (check_bounds && n_far > 0L) {
    stop(sprintf(
      "%s has %d values far outside [%g, %g] (beyond tol = %g).",
      name, n_far, lower, upper, tol
    ))
  }
  
  if (warn_on_clip) {
    n_oob <- sum(below | above)
    warning(sprintf(
      "%s has %d values slightly outside [%g, %g]; clipping applied.",
      name, n_oob, lower, upper
    ))
  }
  
  pmin(pmax(x, lower), upper)
}
