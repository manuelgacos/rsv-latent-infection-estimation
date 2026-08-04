#' Per-subject log-likelihood dispatcher
#'
#' Routes a single subject's contribution to the appropriate likelihood branch
#' (visit vs.\ no-visit) based on \code{subject_i$visit_age}, using the
#' precomputed global and subject-level quantities.
#'
#' This function:
#' \enumerate{
#'   \item Determines whether subject i belongs to I1 (had a visit before age 1)
#'         or I2 (no visit before age 1) from \code{subject_i$visit_age}.
#'   \item Calls the corresponding branch function:
#'         \code{loglik_i_visit()} or \code{loglik_i_no_visit()},
#'         which return the \emph{unweighted} log-likelihood contribution.
#'   \item Multiplies that contribution by \code{subject_i$weight} and returns
#'         the weighted log-likelihood for subject i.
#' }
#'
#' @param subject_i A single subject object created by \code{make_subject()},
#'   containing at least fields \code{visit_age} and \code{weight}.
#'
#' @param subj_pre A list of subject-level precomputations for this subject,
#'   typically the output of \code{rsv_precompute_subject(subject_i, model, control)}.
#'   Must contain at least:
#'   \itemize{
#'     \item \code{lambda_i} – numeric vector of length 365 with
#'           \eqn{\lambda(m + B_i)}, and
#'     \item \code{V_i}      – numeric J x 366 matrix with columns
#'           \eqn{v_i(0),\dots,v_i(365)}.
#'   }
#'
#' @param glob A list of global precomputations for the current parameter values
#'   \code{beta}, \code{eta}, typically the output of
#'   \code{rsv_precompute_global(beta, eta, model, control)}. Must contain:
#'   \itemize{
#'     \item \code{w_vec} – numeric length-365 vector with \eqn{w(m; \beta)},
#'     \item \code{c_vec} – numeric length-365 vector with \eqn{c(m; \eta)}.
#'   }
#'
#' @param model An \code{"rsv_model"} object created by \code{make_rsv_model()}.
#'   Passed through to the branch functions for access to components such as
#'   \code{S_day}.
#'
#' @param beta Numeric vector of length J: parameter vector for the infection
#'   age-density \eqn{w(a; \beta)}. Needed by the branch functions (e.g., via
#'   \code{H_i(beta, V_i)} or \code{Fbar_i(beta, V_i)}).
#'
#' @param eta Numeric vector of length K: parameter vector for the detection
#'   function \eqn{c(a; \eta)}. Needed by the no-visit branch (via
#'   \code{tau_i(U_i, eta)}), and optionally by the visit branch if you choose
#'   to incorporate normalization terms involving \eqn{g_i(\beta,\eta)}.
#'
#' @param control An \code{"rsv_control"} object created by
#'   \code{make_rsv_control()}, passed through to the branch functions to govern
#'   numerical clipping and bound checking.
#'
#' @return Numeric scalar: the \emph{weighted} log-likelihood contribution of
#'   subject i for the current parameter values \code{beta}, \code{eta}.
#'
#' @details
#' This function does not perform any heavy numerical work itself; it simply
#' inspects \code{subject_i$visit_age} and delegates to either
#' \code{loglik_i_visit()} (I1: visit before age 1) or
#' \code{loglik_i_no_visit()} (I2: no visit before age 1). The branch functions
#' are expected to return unweighted log-likelihoods, which are then multiplied
#' by \code{subject_i$weight} here.
#'
#' @export
loglik_i_dispatch <- function(subject_i,
                              subj_pre,
                              glob,
                              model,
                              beta,
                              eta,
                              control) {
  # ---- Cheap structural checks (avoid heavy work inside optimizer) ----------
  
  # visit_age + weight from the subject object
  visit_age <- subject_i$visit_age
  weight    <- subject_i$weight
  
  if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight)) {
    stop("subject_i$weight must be a single finite numeric value.")
  }
  
  # Basic sanity on global precompute
  if (is.null(glob$w_vec) || is.null(glob$c_vec)) {
    stop("`glob` must contain components `w_vec` and `c_vec` from rsv_precompute_global().")
  }
  w_vec <- glob$w_vec
  c_vec <- glob$c_vec
  
  if (!is.numeric(w_vec) || length(w_vec) != 365L ||
      !is.numeric(c_vec) || length(c_vec) != 365L) {
    stop("glob$w_vec and glob$c_vec must be numeric vectors of length 365.")
  }
  
  # Basic sanity on subject precompute
  if (is.null(subj_pre$lambda_i) || is.null(subj_pre$V_i)) {
    stop("`subj_pre` must contain `lambda_i` and `V_i` from rsv_precompute_subject().")
  }
  lambda_i <- subj_pre$lambda_i
  V_i      <- subj_pre$V_i
  
  # ---- Dispatch to visit or no-visit branch ---------------------------------
  
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
    # I1: visit before age 1; enforce that visit_age is in the valid grid
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
  
  # ---- Apply subject weight and return --------------------------------------
  
  if (!is.numeric(ell_core) || length(ell_core) != 1L || !is.finite(ell_core)) {
    stop("Branch function must return a single finite numeric log-likelihood.")
  }
  
  weight * ell_core
}


#' Per-subject log-likelihood: no-visit case (I2)
#'
#' Computes the unweighted log-likelihood contribution for a subject i who
#' had no visit in the first year of life (i.e., \code{visit_age = NA}).
#'
#' In the discrete-day approximation, this contribution is
#' \deqn{
#'   \ell_i(\beta,\eta)
#'     = \log \tau_i(\beta,\eta),
#'   \quad
#'   \tau_i(\beta,\eta) = 1 - \eta^\top U_i(\beta),
#' }
#' where
#' \deqn{
#'   U_i(\beta) = \sum_{m=1}^{365} s(m)\,Q_i(m;\beta),
#'   \quad
#'   Q_i(m;\beta) = \bar F_i(m-1;\beta)\,\pi_i(m;\beta).
#' }
#'
#' This function implements the chain
#' \code{Fbar_i -> pi_from_w -> Q_i -> U_i -> tau_i} using the helpers
#' defined elsewhere in this package.
#'
#' @param lambda_i Numeric vector of length 365 with subject-shifted circulation
#'   values \eqn{\lambda_i[m] = \lambda(m + B_i)} for \eqn{m = 1,\dots,365},
#'   typically obtained from \code{rsv_precompute_subject()}.
#'
#' @param V_i Numeric J x 366 matrix of cumulative kernels for this subject,
#'   as returned by \code{v_i(lambda_i, phi)} inside
#'   \code{rsv_precompute_subject()}. Column 1 corresponds to \eqn{v_i(0)},
#'   columns 2..366 to \eqn{v_i(1)}, \dots, \eqn{v_i(365)}.
#'
#' @param w_vec Numeric vector of length 365 with \eqn{w(m;\beta)} for
#'   \eqn{m = 1,\dots,365}, typically from \code{rsv_precompute_global()}.
#'
#' @param c_vec Numeric vector of length 365 with \eqn{c(m;\eta)}. It is
#'   included for symmetry with the visit-branch interface but is not used
#'   directly in this no-visit branch.
#'
#' @param beta Numeric parameter vector of length J for the infection-age
#'   density \eqn{w(a;\beta)}. Used here via \code{Fbar_i(beta, V_i, ...)}.
#'
#' @param eta Numeric parameter vector of length K for the detection function
#'   \eqn{c(a;\eta)}. Used in \code{tau_i(U_i, eta, ...)}.
#'
#' @param model An \code{"rsv_model"} object created by \code{make_rsv_model()}.
#'   Must contain a numeric matrix \code{S_day} of size K x 365 whose columns
#'   are the day-basis vectors \eqn{s(m)}.
#'
#' @param control An \code{"rsv_control"} object created by
#'   \code{make_rsv_control()}, whose fields \code{check_bounds},
#'   \code{warn_on_clip}, and \code{eps_tau} are propagated to internal
#'   helper functions for numerical clipping and diagnostics.
#'
#' @return A single numeric scalar: the unweighted log-likelihood contribution
#'   \eqn{\ell_i(\beta,\eta) = \log \tau_i(\beta,\eta)} for this subject under
#'   the no-visit case (I2).
#'
#' @export
loglik_i_no_visit <- function(lambda_i,
                              V_i,
                              w_vec,
                              c_vec,
                              beta,
                              eta,
                              model,
                              control) {
  # ---- Basic shape / type checks (cheap) ------------------------------------
  
  # lambda_i
  if (!is.numeric(lambda_i) || length(lambda_i) != 365L) {
    stop("lambda_i must be a numeric vector of length 365.")
  }
  
  # V_i
  if (!is.matrix(V_i) || !is.numeric(V_i)) {
    stop("V_i must be a numeric matrix.")
  }
  if (ncol(V_i) != 366L) {
    stop(sprintf("V_i must have 366 columns (v_i(0..365)); found %d.", ncol(V_i)))
  }
  
  # beta vs V_i rows
  if (!is.numeric(beta) || length(beta) != nrow(V_i)) {
    stop("Length of beta must equal nrow(V_i).")
  }
  
  # w_vec and c_vec
  if (!is.numeric(w_vec) || length(w_vec) != 365L) {
    stop("w_vec must be a numeric vector of length 365.")
  }
  if (!is.numeric(c_vec) || length(c_vec) != 365L) {
    stop("c_vec must be a numeric vector of length 365.")
  }
  # c_vec is not used directly here, but we validate for consistency.
  
  # eta and model$S_day
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
  
  # control
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("control must be an 'rsv_control' object (or compatible list).")
  }
  check_bounds <- isTRUE(control$check_bounds)
  warn_on_clip <- isTRUE(control$warn_on_clip)
  eps_tau      <- if (!is.null(control$eps_tau)) control$eps_tau else 1e-12
  
  # ---- 1) Subject-specific survival: Fbar_i(beta, V_i) ----------------------
  Fbar <- Fbar_i(
    beta        = beta,
    V_i         = V_i,
    include_day0 = TRUE,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  # Fbar should be length 366: d = 0..365
  if (!is.numeric(Fbar) || length(Fbar) != 366L) {
    stop("Fbar_i() must return a numeric vector of length 366 when include_day0 = TRUE.")
  }
  
  # ---- 2) Per-day infection probabilities: pi_from_w ------------------------
  pi <- pi_from_w(
    lambda_shift_i = lambda_i,
    w              = w_vec,
    check_bounds   = check_bounds,
    warn_on_clip   = warn_on_clip
  )
  if (!is.numeric(pi) || length(pi) != 365L) {
    stop("pi_from_w() must return a numeric vector of length 365.")
  }
  
  # ---- 3) Visit kernel Q_i(m; beta) -----------------------------------------
  Q <- Q_i(
    Fbar_i      = Fbar,
    pi_i        = pi,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  if (!is.numeric(Q) || length(Q) != 365L) {
    stop("Q_i() must return a numeric vector of length 365.")
  }
  
  # ---- 4) Aggregated kernel U_i(beta) ---------------------------------------
  U <- U_i(
    Q_i         = Q,
    S_day       = S_day,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  if (!is.numeric(U) || length(U) != length(eta)) {
    stop("U_i() must return a numeric vector of the same length as eta.")
  }
  
  # ---- 5) tau_i(beta, eta) and log-likelihood -------------------------------
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


#' Per-subject log-likelihood: visit case (I1)
#'
#' Computes the unweighted log-likelihood contribution for a subject i who
#' had a bronchiolitis visit in the first year of life (i.e., \code{visit_age}
#' is an integer in 1:365).
#'
#' In the discrete-day approximation, this contribution is
#' \deqn{
#'   \ell_i^{\text{visit}}(\beta,\eta)
#'     = \log c(L_i; \eta)
#'       + \log w(L_i; \beta)
#'       - \beta^\top v_i(L_i),
#' }
#' where:
#' \itemize{
#'   \item \eqn{L_i} is the visit age in days,
#'   \item \eqn{c(L_i; \eta)} is taken from \code{c_vec[L_i]},
#'   \item \eqn{w(L_i; \beta)} is taken from \code{w_vec[L_i]}, and
#'   \item \eqn{v_i(L_i)} is the column \code{V_i[, L_i + 1]} of the cumulative
#'         kernel matrix \code{V_i}.
#' }
#'
#' @section Preconditions on inputs:
#' This function is designed to be used \emph{only} in conjunction with the
#' package's precomputation helpers:
#' \itemize{
#'   \item \code{w_vec} and \code{c_vec} are expected to come from
#'         \code{rsv_precompute_global()}, which enforces nonnegativity via
#'         \code{.enforce_nonneg()} and log-safety via
#'         \code{.stabilize_for_log()}.
#'   \item \code{V_i} is expected to come from \code{rsv_precompute_subject()},
#'         which builds it via \code{v_i(lambda_i, phi)}.
#' }
#' If these preconditions are violated and \code{w_vec[L_i]} or
#' \code{c_vec[L_i]} are nonpositive or non-finite, then:
#' \itemize{
#'   \item if \code{control$check_bounds = TRUE}, an error is thrown;
#'   \item otherwise, a warning is issued (if \code{control$warn_on_clip = TRUE})
#'         and the offending value(s) are repaired to at least
#'         \code{control$eps_log} before taking logs.
#' }
#'
#' @param visit_age Integer in 1:365 giving the first-visit age \eqn{L_i}.
#' @param lambda_i Numeric length-365 vector \eqn{\lambda_i[m] = \lambda(m + B_i)}.
#'   Included for interface symmetry with \code{loglik_i_no_visit()}, but not
#'   used directly in this branch (its effect is already encoded in \code{V_i}).
#' @param V_i Numeric J x 366 matrix of cumulative kernels for this subject,
#'   as returned by \code{rsv_precompute_subject()}. Column 1 corresponds to
#'   \eqn{v_i(0)}, and column d+1 to \eqn{v_i(d)} for d = 1,\dots,365.
#' @param w_vec Numeric length-365 vector with \eqn{w(m; \beta)}, expected
#'   to be the stabilized output of \code{rsv_precompute_global()}.
#' @param c_vec Numeric length-365 vector with \eqn{c(m; \eta)}, expected
#'   to be the stabilized output of \code{rsv_precompute_global()}.
#' @param beta Numeric length-J parameter vector for the infection-age density.
#' @param eta Numeric length-K parameter vector for the detection function.
#'   Included for symmetry with the no-visit branch; currently not used directly.
#' @param model An \code{"rsv_model"} object. Included for interface symmetry
#'   and potential future extensions; not used directly in this branch.
#' @param control An \code{"rsv_control"} object created by
#'   \code{make_rsv_control()}, whose fields \code{check_bounds},
#'   \code{warn_on_clip}, and \code{eps_log} govern diagnostics and repairs
#'   when inconsistencies are detected in \code{w_vec} or \code{c_vec}.
#'
#' @return A single numeric scalar giving the unweighted log-likelihood
#'   contribution \eqn{\ell_i^{\text{visit}}(\beta,\eta)} for this subject.
#'
#' @export
loglik_i_visit <- function(visit_age,
                           lambda_i,
                           V_i,
                           w_vec,
                           c_vec,
                           beta,
                           eta,
                           model,
                           control) {
  # ---- Basic shape / type checks --------------------------------------------
  
  # visit_age
  if (!is.numeric(visit_age) || length(visit_age) != 1L ||
      !is.finite(visit_age) ||
      visit_age < 1L || visit_age > 365L) {
    stop("visit_age must be a single numeric value in 1:365.")
  }
  visit_age <- as.integer(visit_age)
  
  # V_i and beta
  if (!is.matrix(V_i) || !is.numeric(V_i)) {
    stop("V_i must be a numeric matrix.")
  }
  if (ncol(V_i) != 366L) {
    stop(sprintf("V_i must have 366 columns (v_i(0..365)); found %d.", ncol(V_i)))
  }
  if (!is.numeric(beta) || length(beta) != nrow(V_i)) {
    stop("Length of beta must equal nrow(V_i).")
  }
  
  # w_vec and c_vec
  if (!is.numeric(w_vec) || length(w_vec) != 365L) {
    stop("w_vec must be a numeric vector of length 365.")
  }
  if (!is.numeric(c_vec) || length(c_vec) != 365L) {
    stop("c_vec must be a numeric vector of length 365.")
  }
  
  # control object and flags
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("control must be an 'rsv_control' object (or compatible list).")
  }
  check_bounds <- isTRUE(control$check_bounds)
  warn_on_clip <- isTRUE(control$warn_on_clip)
  eps_log      <- if (!is.null(control$eps_log)) control$eps_log else 1e-12
  
  # ---- Extract day-specific weights -----------------------------------------
  
  w_L <- w_vec[visit_age]
  c_L <- c_vec[visit_age]
  
  # Check for nonpositive or non-finite values: these should not occur if
  # w_vec / c_vec came from rsv_precompute_global(), but we guard anyway.
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
  
  # ---- Integrated hazard up to L_i: beta^T v_i(L_i) -------------------------
  
  v_L <- V_i[, visit_age + 1L]
  if (!is.numeric(v_L) || length(v_L) != length(beta)) {
    stop("v_i(L_i) extraction failed: V_i[, visit_age + 1] must be numeric length J.")
  }
  
  H_L <- sum(beta * v_L)
  if (!is.finite(H_L)) {
    stop("Non-finite integrated hazard H_L in loglik_i_visit(); check inputs.")
  }
  
  # ---- Assemble log-likelihood contribution ---------------------------------
  
  ell <- log(c_L) + log(w_L) - H_L
  
  if (!is.numeric(ell) || length(ell) != 1L || !is.finite(ell)) {
    stop("loglik_i_visit() produced a non-finite log-likelihood; check inputs.")
  }
  
  ell
}


#' Total log-likelihood for the RSV infection-age model
#'
#' Computes the total (weighted) log-likelihood
#' \deqn{
#'   \ell(\beta, \eta) = \sum_{i} w_i \, \ell_i(\beta, \eta),
#' }
#' where each subject contribution \eqn{\ell_i} is obtained via
#' \code{loglik_i_dispatch()} and the weights \eqn{w_i} are taken from
#' \code{subject_i$weight}.
#'
#' This function is the main entry point for evaluating the likelihood at a
#' given parameter pair \code{(beta, eta)}. It:
#' \enumerate{
#'   \item Performs global precomputations that depend only on \code{beta},
#'         \code{eta}, and the model (via \code{rsv_precompute_global()}).
#'   \item Ensures that subject-level precomputations are available, either
#'         by using a user-supplied \code{subj_pre_list} or by computing them
#'         on the fly via \code{rsv_precompute_subject()}.
#'   \item Loops over all subjects in \code{data$subjects}, calling
#'         \code{loglik_i_dispatch()} to obtain each subject's (weighted)
#'         log-likelihood contribution.
#'   \item Returns either the total log-likelihood or, optionally, both the
#'         total and the per-subject contributions.
#' }
#'
#' @param beta Numeric vector of length J: parameter vector for the infection
#'   age-density \eqn{w(a; \beta)}.
#'
#' @param eta Numeric vector of length K: parameter vector for the detection
#'   function \eqn{c(a; \eta)}.
#'
#' @param data An \code{"rsv_data"} object created by \code{make_rsv_data()},
#'   containing at least a list component \code{subjects} with one entry per
#'   subject.
#'
#' @param model An \code{"rsv_model"} object created by \code{make_rsv_model()},
#'   providing the day-level basis matrices and circulation curve needed by the
#'   likelihood (e.g., \code{B_day}, \code{S_day}, \code{lambda_global},
#'   \code{phi_day}).
#'
#' @param control An \code{"rsv_control"} object created by
#'   \code{make_rsv_control()}, controlling numerical clipping, bound checks,
#'   and diagnostic behavior. If omitted, a default control object is created.
#'
#' @param subj_pre_list Optional list of subject-level precomputations, typically
#'   created via:
#'   \preformatted{
#'     subj_pre_list <- lapply(data$subjects, rsv_precompute_subject,
#'                             model = model, control = control)
#'   }
#'   Each element should be a list containing at least \code{lambda_i} and
#'   \code{V_i} for the corresponding subject. If \code{subj_pre_list} is
#'   \code{NULL}, subject-level precomputations are performed on the fly.
#'
#' @param return_by_subject Logical; if \code{FALSE} (default), the function
#'   returns a single numeric scalar equal to the total log-likelihood. If
#'   \code{TRUE}, the function returns a list with components:
#'   \describe{
#'     \item{\code{total}}{Total (weighted) log-likelihood.}
#'     \item{\code{per_subject}}{Numeric vector of per-subject (weighted)
#'           log-likelihood contributions, in the same order as
#'           \code{data$subjects}.}
#'   }
#'
#' @return If \code{return_by_subject = FALSE}, a single numeric scalar giving
#'   the total (weighted) log-likelihood. If \code{return_by_subject = TRUE}, a
#'   list with components \code{total} and \code{per_subject} as described
#'   above.
#'
#' @export
rsv_loglik <- function(beta,
                       eta,
                       data,
                       model,
                       control = make_rsv_control(),
                       subj_pre_list = NULL,
                       return_by_subject = FALSE) {
  # ---- Basic argument checks (cheap) ----------------------------------------
  
  # beta, eta
  if (!is.numeric(beta) || any(!is.finite(beta))) {
    stop("rsv_loglik(): 'beta' must be a numeric vector with all finite entries.")
  }
  if (!is.numeric(eta) || any(!is.finite(eta))) {
    stop("rsv_loglik(): 'eta' must be a numeric vector with all finite entries.")
  }
  
  # data
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
  
  # model
  if (is.null(model)) {
    stop("rsv_loglik(): 'model' must be provided (an 'rsv_model' object).")
  }
  
  # control: allow either 'rsv_control' or compatible list
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("rsv_loglik(): 'control' must be an 'rsv_control' object (or compatible list).")
  }
  
  # ---- Global precomputations (depend on beta, eta, model) ------------------
  
  glob <- rsv_precompute_global(
    beta    = beta,
    eta     = eta,
    model   = model,
    control = control
  )
  
  # Expect at least w_vec and c_vec; rsv_precompute_global() should enforce this,
  # but we add a light check for clearer error messages.
  if (is.null(glob$w_vec) || is.null(glob$c_vec)) {
    stop("rsv_loglik(): rsv_precompute_global() did not return 'w_vec' and 'c_vec'.")
  }
  
  # ---- Subject-level precomputations (lambda_i, V_i per subject) ------------
  
  if (is.null(subj_pre_list)) {
    # Compute on the fly
    subj_pre_list <- lapply(
      X   = subjects,
      FUN = rsv_precompute_subject,
      model   = model,
      control = control
    )
  } else {
    # Validate supplied list
    if (!is.list(subj_pre_list) || length(subj_pre_list) != n_subj) {
      stop("rsv_loglik(): 'subj_pre_list' must be a list of length length(data$subjects).")
    }
  }
  
  # ---- Loop over subjects: delegate to dispatcher ---------------------------
  
  per_subject <- numeric(n_subj)
  
  for (i in seq_len(n_subj)) {
    subject_i <- subjects[[i]]
    subj_pre  <- subj_pre_list[[i]]
    
    # loglik_i_dispatch() will perform additional structural checks on subj_pre
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
  
  # ---- Return ----------------------------------------------------------------
  
  if (!return_by_subject) {
    return(ell_total)
  }
  
  # Optionally attach subject IDs as names if they exist consistently
  if (all(vapply(subjects, function(s) !is.null(s$id), logical(1)))) {
    names(per_subject) <- vapply(subjects, function(s) as.character(s$id), character(1))
  }
  
  list(
    total       = ell_total,
    per_subject = per_subject
  )
}


#' Total negative log-likelihood for the RSV infection-age model
#'
#' Computes the negative of the total (weighted) log-likelihood returned by
#' \code{rsv_loglik()}. This is a convenience wrapper intended for use with
#' optimization routines that perform minimization.
#'
#' Formally, if
#' \deqn{
#'   \ell(\beta, \eta) = \sum_i w_i \, \ell_i(\beta, \eta)
#' }
#' is the total (weighted) log-likelihood, then
#' \deqn{
#'   L(\beta, \eta) = -\ell(\beta, \eta)
#' }
#' is the total negative log-likelihood returned by this function when
#' \code{return_by_subject = FALSE}.
#'
#' When \code{return_by_subject = TRUE}, the function returns the negatives of
#' both the total and the per-subject contributions, which can be helpful when
#' inspecting the objective function at the per-subject level in a minimization
#' context.
#'
#' @inheritParams rsv_loglik
#'
#' @param return_by_subject Logical; if \code{FALSE} (default), the function
#'   returns a single numeric scalar equal to the total negative log-likelihood.
#'   If \code{TRUE}, it returns a list with components:
#'   \describe{
#'     \item{\code{total}}{Total negative log-likelihood.}
#'     \item{\code{per_subject}}{Numeric vector of per-subject negative
#'           log-likelihood contributions, in the same order as
#'           \code{data$subjects}.}
#'   }
#'
#' @return If \code{return_by_subject = FALSE}, a single numeric scalar giving
#'   the total negative log-likelihood. If \code{return_by_subject = TRUE}, a
#'   list with components \code{total} and \code{per_subject}, both negated
#'   relative to the output of \code{rsv_loglik()}.
#'
#' @export
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
    # res is a scalar total log-likelihood
    return(-res)
  }
  
  # res is a list(total, per_subject) from rsv_loglik()
  list(
    total       = -res$total,
    per_subject = -res$per_subject
  )
}


#' Objective function: negative log-likelihood with respect to beta
#'
#' Computes the total (weighted) negative log-likelihood
#' \eqn{-\ell(\beta,\eta)} as a scalar function of the infection-age
#' parameter vector \eqn{\beta}. This is a thin wrapper around
#' \code{\link{rsv_negloglik}} intended for use as the \code{fn} argument
#' in \code{\link[stats]{optim}} and related optimization routines.
#'
#' @details
#' This function assumes that all subject-level precomputations
#' (e.g., shifted circulation curves and cumulative kernels) have already
#' been performed and are supplied via \code{subj_pre_list}. No subject-level
#' quantities are recomputed internally.
#'
#' The visit and no-visit likelihood contributions are evaluated using the
#' same numerical safeguards and clipping rules as the full likelihood.
#' The returned value is always a single numeric scalar suitable for
#' gradient-based optimization.
#'
#' @param beta Numeric vector of length \eqn{J}. Infection-age basis coefficients
#'   to be optimized.
#' @param eta Numeric vector of length \eqn{K}. Visit-age basis coefficients,
#'   treated as fixed.
#' @param data An \code{"rsv_data"} object containing the subject list.
#' @param model An \code{"rsv_model"} object containing basis matrices and
#'   circulation curves.
#' @param rsv_control An \code{"rsv_control"} object controlling numerical
#'   checks and clipping behavior. Named \code{rsv_control} to avoid
#'   conflicts with the \code{control} argument of \code{\link[stats]{optim}}.
#' @param subj_pre_list List of subject-level precomputations, typically the
#'   output of \code{\link{rsv_precompute_subject}} for each subject.
#'   Its length must match \code{length(data$subjects)}.
#'
#' @return A numeric scalar giving the total negative log-likelihood
#'   \eqn{-\ell(\beta,\eta)}.
#'
#' @seealso
#'   \code{\link{grad_negloglik_beta}},
#'   \code{\link{rsv_loglik}},
#'   \code{\link{rsv_negloglik}},
#'   \code{\link{rsv_precompute_subject}}
#'
#' @export
obj_negloglik_beta <- function(beta,
                               eta,
                               data,
                               model,
                               rsv_control,
                               subj_pre_list) {
  # --- basic checks (fail fast, cheap) ---
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


#' Objective function: negative log-likelihood with respect to eta
#'
#' Computes the total (weighted) negative log-likelihood
#' \eqn{-\ell(\beta,\eta)} as a scalar function of the visit-age
#' parameter vector \eqn{\eta}, holding \eqn{\beta} fixed.
#'
#' @param eta Numeric vector of length K (parameter to optimize).
#' @param beta Numeric vector of length J (held fixed).
#' @param data rsv_data
#' @param model rsv_model
#' @param rsv_control rsv_control object
#' @param subj_pre_list list of subject-level precomputations
#'
#' @return numeric scalar negative log-likelihood
#' @export
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


#' Quadratic penalty value: (alpha/2) * theta^T M theta
#'
#' @param theta Numeric parameter vector.
#' @param M Symmetric penalty matrix (typically t(D) %*% D).
#' @param alpha Nonnegative smoothing parameter.
#'
#' @return Numeric scalar penalty value.
#' @export
penalty_value <- function(theta, M, alpha) {
  # Fast exit
  if (alpha == 0) {
    return(0.0)
  }
  
  # Minimal checks (cheap, safe)
  stopifnot(is.numeric(theta), is.numeric(alpha), alpha >= 0)
  stopifnot(is.matrix(M), ncol(M) == length(theta), nrow(M) == length(theta))
  
  # Compute 0.5 * alpha * theta^T M theta
  val <- 0.5 * alpha * as.numeric(crossprod(theta, M %*% theta))
  
  val
}


#' Penalized objective: negative log-likelihood with respect to beta
#'
#' Computes the total (weighted) **penalized negative log-likelihood**
#' \eqn{-\ell(\beta,\eta) + \alpha_\beta \, P(\beta)} as a scalar function
#' of the infection-age parameter vector \eqn{\beta}, holding \eqn{\eta} fixed.
#'
#' @details
#' This function augments the negative log-likelihood with a quadratic
#' roughness penalty of the form
#' \deqn{
#'   P(\beta)
#'   =
#'   \frac{1}{2}\beta^\top M_\beta \beta,
#' }
#' where \eqn{M_\beta = D_2^\top D_2} is the second-difference penalty
#' matrix stored in \code{model$M_beta}, and \eqn{\alpha_\beta \ge 0}
#' is the smoothing parameter.
#'
#' The full objective is
#' \deqn{
#'   \tilde L(\beta)
#'   =
#'   -\ell(\beta,\eta)
#'   +
#'   \alpha_\beta \frac{1}{2}\beta^\top M_\beta \beta.
#' }
#'
#' @param beta Numeric vector of length \eqn{J}. Infection-age basis coefficients.
#' @param eta Numeric vector of length \eqn{K}. Visit-age basis coefficients (held fixed).
#' @param alpha_beta Nonnegative scalar smoothing parameter for \eqn{\beta}.
#' @param data rsv_data object.
#' @param model rsv_model object containing \code{M_beta}.
#' @param rsv_control rsv_control object governing numerical safeguards.
#' @param subj_pre_list List of subject-level precomputations.
#'
#' @return Numeric scalar penalized negative log-likelihood.
#'
#' @seealso
#'   \code{\link{obj_grad_pen_negloglik_beta}},
#'   \code{\link{obj_negloglik_beta}},
#'   \code{\link{penalty_value}}
#'
#' @export
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
  
  # Negative log-likelihood
  nll <- obj_negloglik_beta(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  # Penalty
  pen <- penalty_value(
    theta = beta,
    M     = model$M_beta,
    alpha = alpha_beta
  )
  
  nll + pen
}


#' Penalized objective: negative log-likelihood with respect to eta
#'
#' Computes the total (weighted) **penalized negative log-likelihood**
#' \eqn{-\ell(\beta,\eta) + \alpha_\eta \, P(\eta)} as a scalar function
#' of the visit-age parameter vector \eqn{\eta}, holding \eqn{\beta} fixed.
#'
#' @details
#' This function augments the negative log-likelihood with a quadratic
#' roughness penalty of the form
#' \deqn{
#'   P(\eta)
#'   =
#'   \frac{1}{2}\eta^\top M_\eta \eta,
#' }
#' where \eqn{M_\eta = D_2^\top D_2} is the second-difference penalty
#' matrix stored in \code{model$M_eta}, and \eqn{\alpha_\eta \ge 0}
#' is the smoothing parameter.
#'
#' The full objective is
#' \deqn{
#'   \tilde L(\eta)
#'   =
#'   -\ell(\beta,\eta)
#'   +
#'   \alpha_\eta \frac{1}{2}\eta^\top M_\eta \eta.
#' }
#'
#' All likelihood components are evaluated using the same numerical
#' safeguards and clipping rules as \code{\link{rsv_negloglik}}.
#'
#' @param eta Numeric vector of length \eqn{K}. Visit-age basis coefficients
#'   to be optimized.
#'
#' @param beta Numeric vector of length \eqn{J}. Infection-age basis
#'   coefficients, treated as fixed.
#'
#' @param alpha_eta Nonnegative scalar smoothing parameter multiplying
#'   the quadratic penalty.
#'
#' @param data An \code{"rsv_data"} object containing the subject list.
#'
#' @param model An \code{"rsv_model"} object containing basis matrices,
#'   circulation curves, and penalty matrix \code{M_eta}.
#'
#' @param rsv_control An \code{"rsv_control"} object governing numerical
#'   safeguards.
#'
#' @param subj_pre_list List of subject-level precomputations.
#'
#' @return A numeric scalar giving the penalized negative log-likelihood.
#'
#' @seealso
#'   \code{\link{obj_grad_pen_negloglik_eta}},
#'   \code{\link{obj_negloglik_eta}},
#'   \code{\link{penalty_value}},
#'   \code{\link{rsv_negloglik}}
#'
#' @export
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
  
  # Negative log-likelihood
  nll <- obj_negloglik_eta(
    eta           = eta,
    beta          = beta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  # Penalty
  pen <- penalty_value(
    theta = eta,
    M     = model$M_eta,
    alpha = alpha_eta
  )
  
  nll + pen
}


#' Compute the full penalized negative log-likelihood
#'
#' Computes the joint penalized objective function
#' \eqn{
#'   F(\beta, \eta) =
#'   \text{NLL}(\beta, \eta)
#'   + \frac{\alpha}{2}
#'     \left(
#'       \beta^\top M_\beta \beta
#'       +
#'       \eta^\top M_\eta \eta
#'     \right)
#' }
#' where \eqn{\text{NLL}(\beta, \eta)} is the unpenalized negative
#' log-likelihood returned by \code{rsv_negloglik()}.
#'
#' @param beta Numeric vector of spline coefficients for the
#'   infection-age weight function.
#'
#' @param eta Numeric vector of spline coefficients for the
#'   non-birth date covariate effect function.
#'
#' @param alpha Non-negative smoothing parameter controlling
#'   the strength of the quadratic penalties.
#'
#' @param data Data object passed to \code{rsv_negloglik()}.
#'
#' @param model Model object containing penalty matrices
#'   \code{M_beta} and \code{M_eta}.
#'
#' @param rsv_control Control object passed to
#'   \code{rsv_negloglik()}.
#'
#' @param subj_pre_list Precomputed subject-level quantities
#'   used to accelerate likelihood evaluation.
#'
#' @return A numeric scalar giving the value of the full penalized
#'   negative log-likelihood.
#'
#' @seealso \code{\link{rsv_negloglik}},
#'   \code{\link{penalty_value}},
#'   \code{\link{estimate_model_alternating}}
#'
#' @export
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