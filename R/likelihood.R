# Likelihood
#
# Implements subject-level and total likelihood calculations, optimizer
# objectives, and penalized objectives for the RSV model.


#' Dispatch a subject log-likelihood contribution
#'
#' Routes one subject to the visit or no-visit likelihood branch according to
#' `visit_age`, then applies the subject's likelihood weight.
#'
#' @param subject_i A subject record containing `visit_age` and `weight`.
#' @param subj_pre Named list of subject-level precomputations containing
#'   `lambda_i` and `V_i`.
#' @param glob Named list of global precomputations containing `w_vec` and
#'   `c_vec`.
#' @param model An `rsv_model` object passed to the likelihood branch.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param control An `rsv_control` object controlling numerical checks and
#'   stabilization.
#'
#' @return Numeric scalar giving the weighted log-likelihood contribution for
#'   the subject.
loglik_i_dispatch <- function(subject_i,
                              subj_pre,
                              glob,
                              model,
                              beta,
                              eta,
                              control) {
  # Use lightweight checks because this function runs inside optimization.
  visit_age <- subject_i$visit_age
  weight    <- subject_i$weight
  
  if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight)) {
    stop("subject_i$weight must be a single finite numeric value.")
  }
  
  if (is.null(glob$w_vec) || is.null(glob$c_vec)) {
    stop("`glob` must contain components `w_vec` and `c_vec` from rsv_precompute_global().")
  }
  w_vec <- glob$w_vec
  c_vec <- glob$c_vec
  
  if (!is.numeric(w_vec) || length(w_vec) != 365L ||
      !is.numeric(c_vec) || length(c_vec) != 365L) {
    stop("glob$w_vec and glob$c_vec must be numeric vectors of length 365.")
  }
  
  if (is.null(subj_pre$lambda_i) || is.null(subj_pre$V_i)) {
    stop("`subj_pre` must contain `lambda_i` and `V_i` from rsv_precompute_subject().")
  }
  lambda_i <- subj_pre$lambda_i
  V_i      <- subj_pre$V_i
  
  if (is.na(visit_age)) {
    # I2: no visit before age 1
    ell_core <- loglik_i_no_visit(
      lambda_i = lambda_i,
      V_i      = V_i,
      w_vec    = w_vec,
      c_vec    = c_vec,
      beta     = beta,
      eta      = eta,
      model    = model,
      control  = control
    )
    
  } else {
    # I1: visit before age 1; validate the observed visit age.
    if (!is.numeric(visit_age) || length(visit_age) != 1L ||
        !is.finite(visit_age) ||
        visit_age < 1L || visit_age > 365L) {
      stop("subject_i$visit_age must be NA or a single integer in 1:365.")
    }
    visit_age <- as.integer(visit_age)
    
    ell_core <- loglik_i_visit(
      visit_age = visit_age,
      lambda_i  = lambda_i,
      V_i       = V_i,
      w_vec     = w_vec,
      c_vec     = c_vec,
      beta      = beta,
      eta       = eta,
      model     = model,
      control   = control
    )
  }
  
  if (!is.numeric(ell_core) || length(ell_core) != 1L || !is.finite(ell_core)) {
    stop("Branch function must return a single finite numeric log-likelihood.")
  }
  
  weight * ell_core
}


#' Evaluate the no-visit log-likelihood contribution
#'
#' Computes the unweighted log-likelihood contribution for a subject with no
#' healthcare visit during the first year of life.
#'
#' @param lambda_i Numeric vector of length 365 containing the
#'   subject-specific RSV circulation curve over ages 1 through 365 days.
#' @param V_i Numeric `J x 366` matrix containing the cumulative subject-level
#'   kernels. Column 1 represents \eqn{v_i(0)}, and column `d + 1` represents
#'   \eqn{v_i(d)} for days 1 through 365.
#' @param w_vec Numeric vector of length 365 containing the infection-age curve
#'   \eqn{w(a)} over ages 1 through 365 days.
#' @param c_vec Numeric vector of length 365 containing the healthcare-visit
#'   curve \eqn{c(a)}. This argument is validated for interface consistency
#'   but is not used directly in this likelihood branch.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param model An `rsv_model` object containing the healthcare-visit basis
#'   `S_day`.
#' @param control An `rsv_control` object controlling numerical checks and
#'   stabilization.
#'
#' @return Numeric scalar giving the unweighted no-visit log-likelihood
#'   contribution for the subject.
#'
#' @details
#' The no-visit contribution is
#' \deqn{
#'   \ell_i(\beta,\eta)
#'   =
#'   \log \tau_i(\beta,\eta),
#' }
#' where
#' \deqn{
#'   \tau_i(\beta,\eta)
#'   =
#'   1 - \eta^\top U_i(\beta),
#' }
#' with
#' \deqn{
#'   U_i(\beta)
#'   =
#'   \sum_{m=1}^{365} s(m)Q_i(m;\beta),
#'   \qquad
#'   Q_i(m;\beta)
#'   =
#'   \bar F_i(m-1;\beta)\pi_i(m;\beta).
#' }
#'
#' The computation follows
#' `Fbar_i()` -> `pi_from_w()` -> `Q_i()` -> `U_i()` -> `tau_i()`.
loglik_i_no_visit <- function(lambda_i,
                              V_i,
                              w_vec,
                              c_vec,
                              beta,
                              eta,
                              model,
                              control) {
  if (!is.numeric(lambda_i) || length(lambda_i) != 365L) {
    stop("lambda_i must be a numeric vector of length 365.")
  }
  
  if (!is.matrix(V_i) || !is.numeric(V_i)) {
    stop("V_i must be a numeric matrix.")
  }
  if (ncol(V_i) != 366L) {
    stop(sprintf("V_i must have 366 columns (v_i(0..365)); found %d.", ncol(V_i)))
  }
  
  if (!is.numeric(beta) || length(beta) != nrow(V_i)) {
    stop("Length of beta must equal nrow(V_i).")
  }
  
  if (!is.numeric(w_vec) || length(w_vec) != 365L) {
    stop("w_vec must be a numeric vector of length 365.")
  }
  # Validate c_vec for consistency with the visit-branch interface.
  if (!is.numeric(c_vec) || length(c_vec) != 365L) {
    stop("c_vec must be a numeric vector of length 365.")
  }
  
  if (!is.numeric(eta)) {
    stop("eta must be numeric.")
  }
  if (is.null(model$S_day)) {
    stop("model$S_day is missing; it must be provided in the rsv_model object.")
  }
  S_day <- model$S_day
  if (!is.matrix(S_day) && !inherits(S_day, "Matrix")) {
    stop("model$S_day must be a numeric matrix (or Matrix).")
  }
  if (!is.numeric(S_day)) {
    stop("model$S_day must be numeric.")
  }
  if (ncol(S_day) != 365L) {
    stop(sprintf("model$S_day must have 365 columns; found %d.", ncol(S_day)))
  }
  if (length(eta) != nrow(S_day)) {
    stop("Length of eta must equal nrow(model$S_day).")
  }
  
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("control must be an 'rsv_control' object (or compatible list).")
  }
  check_bounds <- isTRUE(control$check_bounds)
  warn_on_clip <- isTRUE(control$warn_on_clip)
  eps_tau      <- if (!is.null(control$eps_tau)) control$eps_tau else 1e-12
  
  # Compute subject-specific survival over days 0 through 365.
  Fbar <- Fbar_i(
    beta        = beta,
    V_i         = V_i,
    include_day0 = TRUE,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  if (!is.numeric(Fbar) || length(Fbar) != 366L) {
    stop("Fbar_i() must return a numeric vector of length 366 when include_day0 = TRUE.")
  }
  
  # Compute the conditional infection probability for each age day.
  pi <- pi_from_w(
    lambda_shift_i = lambda_i,
    w              = w_vec,
    check_bounds   = check_bounds,
    warn_on_clip   = warn_on_clip
  )
  if (!is.numeric(pi) || length(pi) != 365L) {
    stop("pi_from_w() must return a numeric vector of length 365.")
  }
  
  # Compute the infection-day probability mass.
  Q <- Q_i(
    Fbar_i      = Fbar,
    pi_i        = pi,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  if (!is.numeric(Q) || length(Q) != 365L) {
    stop("Q_i() must return a numeric vector of length 365.")
  }
  
  # Aggregate the infection-day mass over the healthcare-visit basis.
  U <- U_i(
    Q_i         = Q,
    S_day       = S_day,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  if (!is.numeric(U) || length(U) != length(eta)) {
    stop("U_i() must return a numeric vector of the same length as eta.")
  }
  
  # Compute the no-visit probability before taking its logarithm.
  tau <- tau_i(
    U_i         = U,
    eta         = eta,
    eps_tau     = eps_tau,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau_i() must return a single positive finite numeric scalar.")
  }
  
  log(tau)
}


#' Evaluate the visit log-likelihood contribution
#'
#' Computes the unweighted log-likelihood contribution for a subject with a
#' healthcare visit during the first year of life.
#'
#' @param visit_age Integer scalar in 1:365 giving the healthcare-visit age
#'   \eqn{L_i}.
#' @param lambda_i Numeric vector of length 365 containing the
#'   subject-specific RSV circulation curve. This argument is included for
#'   consistency with the no-visit likelihood branch but is not used directly.
#' @param V_i Numeric `J x 366` matrix containing the cumulative subject-level
#'   kernels. Column 1 represents \eqn{v_i(0)}, and column `d + 1` represents
#'   \eqn{v_i(d)} for days 1 through 365.
#' @param w_vec Numeric vector of length 365 containing the infection-age curve
#'   \eqn{w(a)} over ages 1 through 365 days.
#' @param c_vec Numeric vector of length 365 containing the healthcare-visit
#'   curve \eqn{c(a)} over ages 1 through 365 days.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients. This argument is included for interface consistency
#'   but is not used directly.
#' @param model An `rsv_model` object included for interface consistency but
#'   not used directly in this likelihood branch.
#' @param control An `rsv_control` object controlling numerical checks and
#'   stabilization.
#'
#' @return Numeric scalar giving the unweighted visit log-likelihood
#'   contribution for the subject.
#'
#' @details
#' For a visit at age \eqn{L_i}, the contribution is
#' \deqn{
#'   \ell_i(\beta,\eta)
#'   =
#'   \log c(L_i;\eta)
#'   +
#'   \log w(L_i;\beta)
#'   -
#'   \beta^\top v_i(L_i).
#' }
#' The cumulative kernel \eqn{v_i(L_i)} is stored in column `L_i + 1` of
#' `V_i` because column 1 represents day zero.
#'
#' If the selected values of `w_vec` or `c_vec` are nonpositive or
#' non-finite, `control` determines whether the function stops or applies
#' log-scale stabilization.
loglik_i_visit <- function(visit_age,
                           lambda_i,
                           V_i,
                           w_vec,
                           c_vec,
                           beta,
                           eta,
                           model,
                           control) {
  if (!is.numeric(visit_age) || length(visit_age) != 1L ||
      !is.finite(visit_age) ||
      visit_age < 1L || visit_age > 365L) {
    stop("visit_age must be a single numeric value in 1:365.")
  }
  visit_age <- as.integer(visit_age)
  
  if (!is.matrix(V_i) || !is.numeric(V_i)) {
    stop("V_i must be a numeric matrix.")
  }
  if (ncol(V_i) != 366L) {
    stop(sprintf("V_i must have 366 columns (v_i(0..365)); found %d.", ncol(V_i)))
  }
  if (!is.numeric(beta) || length(beta) != nrow(V_i)) {
    stop("Length of beta must equal nrow(V_i).")
  }
  
  if (!is.numeric(w_vec) || length(w_vec) != 365L) {
    stop("w_vec must be a numeric vector of length 365.")
  }
  if (!is.numeric(c_vec) || length(c_vec) != 365L) {
    stop("c_vec must be a numeric vector of length 365.")
  }
  
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("control must be an 'rsv_control' object (or compatible list).")
  }
  check_bounds <- isTRUE(control$check_bounds)
  warn_on_clip <- isTRUE(control$warn_on_clip)
  eps_log      <- if (!is.null(control$eps_log)) control$eps_log else 1e-12
  
  w_L <- w_vec[visit_age]
  c_L <- c_vec[visit_age]
  
  # Guard against invalid curve values before taking logarithms.
  bad_w <- !is.finite(w_L) || w_L <= 0
  bad_c <- !is.finite(c_L) || c_L <= 0
  
  if (bad_w || bad_c) {
    msg <- sprintf(
      paste0(
        "Nonpositive or non-finite w_vec[%d] or c_vec[%d] detected in loglik_i_visit(). ",
        "This suggests that rsv_precompute_global() was not used, or that ",
        "its invariants were violated."
      ),
      visit_age, visit_age
    )
    
    if (check_bounds) {
      stop(msg)
    } else {
      if (warn_on_clip) {
        warning(paste0(msg, " Repairing to eps_log for log-safety."))
      }
      if (bad_w) w_L <- eps_log
      if (bad_c) c_L <- eps_log
    }
  }
  
  # Column visit_age + 1 represents v_i(L_i) because column 1 is day zero.
  v_L <- V_i[, visit_age + 1L]
  if (!is.numeric(v_L) || length(v_L) != length(beta)) {
    stop("v_i(L_i) extraction failed: V_i[, visit_age + 1] must be numeric length J.")
  }
  
  H_L <- sum(beta * v_L)
  if (!is.finite(H_L)) {
    stop("Non-finite integrated hazard H_L in loglik_i_visit(); check inputs.")
  }
  
  ell <- log(c_L) + log(w_L) - H_L
  
  if (!is.numeric(ell) || length(ell) != 1L || !is.finite(ell)) {
    stop("loglik_i_visit() produced a non-finite log-likelihood; check inputs.")
  }
  
  ell
}


#' Evaluate the total RSV log-likelihood
#'
#' Computes the total weighted log-likelihood across all subjects for the
#' current infection-age and healthcare-visit spline coefficients.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases and
#'   seasonal RSV circulation curve.
#' @param control An `rsv_control` object controlling numerical checks and
#'   stabilization.
#' @param subj_pre_list Optional list of subject-level precomputations, with
#'   one element per subject. Each element must contain `lambda_i` and `V_i`.
#'   If `NULL`, the precomputations are constructed internally.
#' @param return_by_subject Logical scalar indicating whether to return the
#'   individual weighted subject contributions in addition to their total.
#'
#' @return If `return_by_subject = FALSE`, a numeric scalar giving the total
#'   weighted log-likelihood. If `TRUE`, a named list containing:
#' \describe{
#'   \item{\code{total}}{Numeric scalar giving the total weighted
#'     log-likelihood.}
#'   \item{\code{per_subject}}{Numeric vector containing the weighted
#'     log-likelihood contribution for each subject.}
#' }
#'
#' @details
#' The total log-likelihood is
#' \deqn{
#'   \ell(\beta,\eta)
#'   =
#'   \sum_{i=1}^{n} \omega_i \ell_i(\beta,\eta),
#' }
#' where \eqn{\omega_i} is the subject weight and \eqn{\ell_i} is the
#' appropriate visit or no-visit contribution.
#'
#' Global infection-age and healthcare-visit curves are computed once for the
#' current coefficients. Subject-level precomputations are then reused when
#' evaluating each contribution.
rsv_loglik <- function(beta,
                       eta,
                       data,
                       model,
                       control = make_rsv_control(),
                       subj_pre_list = NULL,
                       return_by_subject = FALSE) {
  
  if (!is.numeric(beta) || any(!is.finite(beta))) {
    stop("rsv_loglik(): 'beta' must be a numeric vector with all finite entries.")
  }
  if (!is.numeric(eta) || any(!is.finite(eta))) {
    stop("rsv_loglik(): 'eta' must be a numeric vector with all finite entries.")
  }
  
  if (is.null(data) || is.null(data$subjects)) {
    stop("rsv_loglik(): 'data' must be an 'rsv_data' object with a 'subjects' component.")
  }
  subjects <- data$subjects
  if (!is.list(subjects)) {
    stop("rsv_loglik(): data$subjects must be a list of subject objects.")
  }
  n_subj <- length(subjects)
  if (n_subj == 0L) {
    stop("rsv_loglik(): data contains no subjects; cannot compute likelihood.")
  }
  
  if (is.null(model)) {
    stop("rsv_loglik(): 'model' must be provided (an 'rsv_model' object).")
  }
  
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("rsv_loglik(): 'control' must be an 'rsv_control' object (or compatible list).")
  }
  
  # Compute parameter-dependent curves once for all subjects.
  glob <- rsv_precompute_global(
    beta    = beta,
    eta     = eta,
    model   = model,
    control = control
  )
  
  # Retain a light check for clearer errors from the precomputation.
  if (is.null(glob$w_vec) || is.null(glob$c_vec)) {
    stop("rsv_loglik(): rsv_precompute_global() did not return 'w_vec' and 'c_vec'.")
  }
  
  # Reuse supplied subject-level precomputations when available.
  if (is.null(subj_pre_list)) {
    subj_pre_list <- lapply(
      X   = subjects,
      FUN = rsv_precompute_subject,
      model   = model,
      control = control
    )
  } else {
    if (!is.list(subj_pre_list) || length(subj_pre_list) != n_subj) {
      stop("rsv_loglik(): 'subj_pre_list' must be a list of length length(data$subjects).")
    }
  }
  
  # Evaluate the appropriate likelihood branch for each subject.
  per_subject <- numeric(n_subj)
  
  for (i in seq_len(n_subj)) {
    subject_i <- subjects[[i]]
    subj_pre  <- subj_pre_list[[i]]
    
    ell_i <- loglik_i_dispatch(
      subject_i = subject_i,
      subj_pre  = subj_pre,
      glob      = glob,
      model     = model,
      beta      = beta,
      eta       = eta,
      control   = control
    )
    
    per_subject[i] <- ell_i
  }
  
  ell_total <- sum(per_subject)
  
  if (!return_by_subject) {
    return(ell_total)
  }
  
  # Use subject IDs as names only when every record provides one.
  if (all(vapply(subjects, function(s) !is.null(s$id), logical(1)))) {
    names(per_subject) <- vapply(subjects, function(s) as.character(s$id), character(1))
  }
  
  list(
    total       = ell_total,
    per_subject = per_subject
  )
}


#' Evaluate the total negative log-likelihood
#'
#' Returns the negative of the weighted log-likelihood computed by
#' `rsv_loglik()`. This form is used by optimization routines that minimize
#' their objective function.
#'
#' @inheritParams rsv_loglik
#'
#' @return If `return_by_subject = FALSE`, a numeric scalar giving the total
#'   negative log-likelihood. If `TRUE`, a named list containing:
#' \describe{
#'   \item{\code{total}}{Numeric scalar giving the total negative
#'     log-likelihood.}
#'   \item{\code{per_subject}}{Numeric vector containing the negative
#'     log-likelihood contribution for each subject.}
#' }
rsv_negloglik <- function(beta,
                          eta,
                          data,
                          model,
                          control = make_rsv_control(),
                          subj_pre_list = NULL,
                          return_by_subject = FALSE) {
  res <- rsv_loglik(
    beta              = beta,
    eta               = eta,
    data              = data,
    model             = model,
    control           = control,
    subj_pre_list     = subj_pre_list,
    return_by_subject = return_by_subject
  )
  
  if (!return_by_subject) {
    return(-res)
  }
  
  # Negate both components when subject-level contributions are requested.
  list(
    total       = -res$total,
    per_subject = -res$per_subject
  )
}


#' Evaluate the negative log-likelihood for beta optimization
#'
#' Returns the total weighted negative log-likelihood as a function of
#' \eqn{\beta}, with \eqn{\eta} held fixed. This wrapper is used by the
#' blockwise optimization routine for the infection-age coefficients.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients to optimize.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients held fixed during optimization.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the model components used
#'   in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks
#'   and stabilization. The argument name avoids conflict with an optimizer's
#'   `control` argument.
#' @param subj_pre_list List of subject-level precomputations, with one
#'   element per subject.
#'
#' @return Numeric scalar giving the total weighted negative log-likelihood.
obj_negloglik_beta <- function(beta,
                               eta,
                               data,
                               model,
                               rsv_control,
                               subj_pre_list) {
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  rsv_negloglik(
    beta = beta,
    eta  = eta,
    data = data,
    model = model,
    control = rsv_control,
    subj_pre_list = subj_pre_list,
    return_by_subject = FALSE
  )
}


#' Evaluate the negative log-likelihood for eta optimization
#'
#' Returns the total weighted negative log-likelihood as a function of
#' \eqn{\eta}, with \eqn{\beta} held fixed. This wrapper is used by the
#' blockwise optimization routine for the healthcare-visit coefficients.
#'
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients to optimize.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients held fixed during optimization.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the model components used
#'   in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks
#'   and stabilization. The argument name avoids conflict with an optimizer's
#'   `control` argument.
#' @param subj_pre_list List of subject-level precomputations, with one
#'   element per subject.
#'
#' @return Numeric scalar giving the total weighted negative log-likelihood.
obj_negloglik_eta <- function(eta,
                              beta,
                              data,
                              model,
                              rsv_control,
                              subj_pre_list) {
  stopifnot(is.numeric(eta), all(is.finite(eta)))
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  rsv_negloglik(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    control       = rsv_control,
    subj_pre_list = subj_pre_list,
    return_by_subject = FALSE
  )
}


#' Evaluate a quadratic smoothing penalty
#'
#' Computes the scaled quadratic penalty for a coefficient vector and its
#' penalty matrix.
#'
#' @param theta Numeric vector of coefficients.
#' @param M Numeric square matrix with dimensions matching `length(theta)`,
#'   containing the quadratic penalty matrix.
#' @param alpha Nonnegative numeric scalar controlling the penalty strength.
#'
#' @return Numeric scalar giving the quadratic penalty value.
#'
#' @details
#' The penalty is
#' \deqn{
#'   \frac{\alpha}{2}\theta^\top M\theta.
#' }
#' When `alpha = 0`, the function returns zero without evaluating the
#' quadratic form.
penalty_value <- function(theta, M, alpha) {
  # Return immediately when no penalty is applied.
  if (alpha == 0) {
    return(0.0)
  }
  
  stopifnot(is.numeric(theta), is.numeric(alpha), alpha >= 0)
  stopifnot(is.matrix(M), ncol(M) == length(theta), nrow(M) == length(theta))
  
  val <- 0.5 * alpha * as.numeric(crossprod(theta, M %*% theta))
  
  val
}


#' Evaluate the penalized objective for beta optimization
#'
#' Computes the total weighted negative log-likelihood plus the quadratic
#' smoothing penalty for \eqn{\beta}, with \eqn{\eta} held fixed.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients to optimize.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients held fixed during optimization.
#' @param alpha_beta Nonnegative numeric scalar controlling the smoothing
#'   penalty for \eqn{\beta}.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the penalty matrix `M_beta`
#'   and the components used in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks
#'   and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one
#'   element per subject.
#'
#' @return Numeric scalar giving the penalized negative log-likelihood.
#'
#' @details
#' The objective is
#' \deqn{
#'   -\ell(\beta,\eta)
#'   +
#'   \frac{\alpha_\beta}{2}
#'   \beta^\top M_\beta \beta.
#' }
obj_pen_negloglik_beta <- function(beta,
                                   eta,
                                   alpha_beta,
                                   data,
                                   model,
                                   rsv_control,
                                   subj_pre_list) {
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.numeric(alpha_beta), length(alpha_beta) == 1L, alpha_beta >= 0)
  stopifnot(!is.null(model$M_beta))
  
  # Compute the unpenalized negative log-likelihood.
  nll <- obj_negloglik_beta(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  # Compute the quadratic smoothing penalty for beta.
  pen <- penalty_value(
    theta = beta,
    M     = model$M_beta,
    alpha = alpha_beta
  )
  
  nll + pen
}


#' Evaluate the penalized objective for eta optimization
#'
#' Computes the total weighted negative log-likelihood plus the quadratic
#' smoothing penalty for \eqn{\eta}, with \eqn{\beta} held fixed.
#'
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients to optimize.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients held fixed during optimization.
#' @param alpha_eta Nonnegative numeric scalar controlling the smoothing
#'   penalty for \eqn{\eta}.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the penalty matrix `M_eta`
#'   and the components used in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks
#'   and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one
#'   element per subject.
#'
#' @return Numeric scalar giving the penalized negative log-likelihood.
#'
#' @details
#' The objective is
#' \deqn{
#'   -\ell(\beta,\eta)
#'   +
#'   \frac{\alpha_\eta}{2}
#'   \eta^\top M_\eta \eta.
#' }
obj_pen_negloglik_eta <- function(eta,
                                  beta,
                                  alpha_eta,
                                  data,
                                  model,
                                  rsv_control,
                                  subj_pre_list) {
  stopifnot(is.numeric(eta), all(is.finite(eta)))
  stopifnot(is.numeric(alpha_eta), length(alpha_eta) == 1L, alpha_eta >= 0)
  stopifnot(!is.null(model$M_eta))
  
  # Compute the unpenalized negative log-likelihood.
  nll <- obj_negloglik_eta(
    eta           = eta,
    beta          = beta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  # Compute the quadratic smoothing penalty for eta.
  pen <- penalty_value(
    theta = eta,
    M     = model$M_eta,
    alpha = alpha_eta
  )
  
  nll + pen
}


#' Evaluate the joint penalized objective
#'
#' Computes the total weighted negative log-likelihood plus quadratic
#' smoothing penalties for the infection-age and healthcare-visit spline
#' coefficients.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param alpha Nonnegative numeric scalar controlling the smoothing penalty
#'   for both coefficient vectors.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the penalty matrices `M_beta`
#'   and `M_eta` and the components used in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks
#'   and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one
#'   element per subject.
#'
#' @return Numeric scalar giving the joint penalized objective value.
#'
#' @details
#' The objective is
#' \deqn{
#'   -\ell(\beta,\eta)
#'   +
#'   \frac{\alpha}{2}
#'   \left(
#'     \beta^\top M_\beta \beta
#'     +
#'     \eta^\top M_\eta \eta
#'   \right).
#' }
full_penalized_objective <- function(
    beta,
    eta,
    alpha,
    data,
    model,
    rsv_control,
    subj_pre_list
) {
  nll <- rsv_negloglik(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    control       = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  pen_beta <- penalty_value(
    theta = beta,
    M     = model$M_beta,
    alpha = alpha
  )
  
  pen_eta <- penalty_value(
    theta = eta,
    M     = model$M_eta,
    alpha = alpha
  )
  
  obj <- nll + pen_beta + pen_eta
  
  return(obj)
}
