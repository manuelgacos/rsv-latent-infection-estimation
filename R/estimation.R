#' Gradient of total log-likelihood with respect to beta (fast version)
#'
#' Computes the analytic gradient of the total (weighted) log-likelihood
#' \eqn{\ell(\beta,\eta)} with respect to the infection-age parameter vector
#' \eqn{\beta}, using a vectorized implementation for the no-visit (I2) case.
#'
#' This function is algebraically equivalent to \code{\link{grad_beta}},
#' but replaces the per-day loop in the I2 term with matrix–vector products,
#' resulting in substantially improved performance for large datasets or
#' large basis dimensions.
#'
#' @details
#' The gradient is computed by summing per-subject contributions,
#' weighted by \code{subject_i$weight}, and split into two cases:
#'
#' \strong{Visit case (I1):}
#' If subject \eqn{i} has a visit at age \eqn{L_i},
#' \deqn{
#'   \nabla_\beta \ell_i
#'   =
#'   \frac{b(L_i)}{w(L_i;\beta)} - v_i(L_i),
#' }
#' where \eqn{b(L_i)} is the day-basis vector, \eqn{w(L_i;\beta) = \beta^\top b(L_i)},
#' and \eqn{v_i(L_i)} is the cumulative kernel at age \eqn{L_i}.
#'
#' \strong{No-visit case (I2):}
#' If subject \eqn{i} has no visit in the first year of life,
#' \deqn{
#'   \nabla_\beta \ell_i
#'   =
#'   -\frac{1}{\tau_i(\beta,\eta)}
#'   \sum_{m=1}^{365}
#'   c(m;\eta)
#'   \left[
#'     \bar F_i(m;\beta)\,\lambda_i(m)\,b(m)
#'     -
#'     Q_i(m;\beta)\,v_i(m-1)
#'   \right],
#' }
#' where:
#' \itemize{
#'   \item \eqn{\bar F_i(m;\beta)} is the survival probability up to day \eqn{m},
#'   \item \eqn{\lambda_i(m)} is the subject-specific circulation curve,
#'   \item \eqn{Q_i(m;\beta) = \bar F_i(m-1;\beta)\,\pi_i(m;\beta)} is the
#'         infection-day mass,
#'   \item \eqn{\tau_i(\beta,\eta) = 1 - \sum_m c(m;\eta)\,Q_i(m;\beta)} is the
#'         probability of no visit.
#' }
#'
#' The inner sum is evaluated efficiently using matrix–vector products:
#' \itemize{
#'   \item \code{B_day \%*\% a}, with \eqn{a(m) = c(m)\,\bar F_i(m)\,\lambda_i(m)},
#'   \item \code{V_i[,1:365] \%*\% b}, with \eqn{b(m) = c(m)\,Q_i(m)}.
#' }
#'
#' @section Numerical considerations:
#' This function reuses the same numerical safeguards as the likelihood:
#' \itemize{
#'   \item The vectors \code{w_vec} and \code{c_vec} are obtained from
#'         \code{\link{rsv_precompute_global}}, which enforces nonnegativity
#'         and log-safety via \code{control$eps_log}.
#'   \item The no-visit probability \eqn{\tau_i} is clipped below by
#'         \code{control$eps_tau} to avoid numerical instability.
#'   \item Clipping operations are treated as fixed for differentiation;
#'         the returned gradient corresponds to the analytic derivative of
#'         the unclipped likelihood.
#' }
#'
#' @param beta Numeric vector of length \eqn{J}. Infection-age basis coefficients.
#' @param eta Numeric vector of length \eqn{K}. Visit-age basis coefficients
#'   (held fixed in this gradient).
#' @param data An \code{"rsv_data"} object containing the subject list.
#' @param model An \code{"rsv_model"} object containing the basis matrices
#'   and circulation curve.
#' @param control An \code{"rsv_control"} object controlling numerical checks,
#'   clipping behavior, and tolerance parameters.
#' @param subj_pre_list Optional list of subject-level precomputations,
#'   typically the output of \code{\link{rsv_precompute_subject}} for each subject.
#'   If \code{NULL}, subject-level quantities are computed internally.
#'
#' @return A numeric vector of length \eqn{J} giving the gradient
#'   \eqn{\nabla_\beta \ell(\beta,\eta)} of the total (weighted) log-likelihood.
#'
#' @seealso
#'   \code{\link{grad_beta}},
#'   \code{\link{grad_negloglik_beta}},
#'   \code{\link{rsv_loglik}},
#'   \code{\link{rsv_precompute_global}},
#'   \code{\link{rsv_precompute_subject}}
#'
#' @export
grad_beta_fast <- function(beta,
                           eta,
                           data,
                           model,
                           control = make_rsv_control(),
                           subj_pre_list = NULL) {
  # --- basic checks ---
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.numeric(eta),  all(is.finite(eta)))
  stopifnot(!is.null(data$subjects), is.list(data$subjects))
  stopifnot(inherits(model, "rsv_model"))
  
  check_bounds <- isTRUE(control$check_bounds)
  warn_on_clip <- isTRUE(control$warn_on_clip)
  eps_log      <- if (!is.null(control$eps_log)) control$eps_log else 1e-12
  eps_tau      <- if (!is.null(control$eps_tau)) control$eps_tau else 1e-12
  
  subjects <- data$subjects
  n_subj   <- length(subjects)
  J        <- length(beta)
  
  # --- global precompute: w_vec, c_vec ---
  glob  <- rsv_precompute_global(beta = beta, eta = eta, model = model, control = control)
  w_vec <- glob$w_vec
  c_vec <- glob$c_vec
  
  # --- subject precompute if needed ---
  if (is.null(subj_pre_list)) {
    subj_pre_list <- lapply(subjects, rsv_precompute_subject, model = model, control = control)
  } else {
    stopifnot(is.list(subj_pre_list), length(subj_pre_list) == n_subj)
  }
  
  grad  <- rep(0.0, J)
  B_day <- model$B_day
  
  for (i in seq_len(n_subj)) {
    subj   <- subjects[[i]]
    pre    <- subj_pre_list[[i]]
    weight <- subj$weight
    if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight)) {
      stop("subject weight must be finite scalar.")
    }
    
    lambda_i <- pre$lambda_i        # length 365
    V_i      <- pre$V_i             # J x 366
    
    L_i <- subj$visit_age
    
    if (!is.na(L_i)) {
      # -------------------------
      # I1 term: b(L)/w(L) - v(L)
      # -------------------------
      L_i <- as.integer(L_i)
      b_L <- B_day[, L_i]
      w_L <- w_vec[L_i]
      if (!is.finite(w_L) || w_L <= 0) {
        if (check_bounds) stop("w_vec[L_i] nonpositive/nonfinite in grad_beta I1.")
        if (warn_on_clip) warning("w_vec[L_i] nonpositive/nonfinite in grad_beta I1; repairing.")
        w_L <- eps_log
      }
      v_L <- V_i[, L_i + 1L]
      grad <- grad + weight * ((b_L / w_L) - v_L)
      
    } else {
      # -------------------------
      # I2 term (vectorized inner)
      # -------------------------
      Fbar <- Fbar_i(
        beta         = beta,
        V_i          = V_i,
        include_day0 = TRUE,
        check_bounds = check_bounds,
        warn_on_clip = warn_on_clip
      )
      
      pi <- pi_from_w(
        lambda_shift_i = lambda_i,
        w              = w_vec,
        check_bounds   = check_bounds,
        warn_on_clip   = warn_on_clip
      )
      
      Q <- Q_i(
        Fbar_i       = Fbar,
        pi_i         = pi,
        check_bounds = check_bounds,
        warn_on_clip = warn_on_clip
      )
      
      g_i <- sum(c_vec * Q)
      tau <- 1 - g_i
      tau <- .clip_to_range(
        tau, lower = eps_tau, upper = 1.0, name = "tau_i",
        check_bounds = check_bounds, warn_on_clip = warn_on_clip
      )
      
      # a[m] = c(m) * Fbar(m) * lambda_i(m)  where Fbar(m) = Fbar[m+1]
      a <- c_vec * (Fbar[2:366] * lambda_i)  # length 365
      
      # b[m] = c(m) * Q(m)
      b <- c_vec * Q                         # length 365
      
      # inner = B_day %*% a - V_i[,1:365] %*% b   (both yield length J)
      inner <- as.vector(B_day %*% a) - as.vector(V_i[, 1:365, drop = FALSE] %*% b)
      
      grad <- grad + weight * (-(1 / tau)) * inner
    }
  }
  
  grad
}


#' Gradient of negative log-likelihood w.r.t. beta
#'
#' Wrapper around grad_beta_fast() that returns the gradient of
#' the *negative* log-likelihood.
#'
#' @param beta numeric vector of length J
#' @param eta  numeric vector of length K (held fixed)
#' @param data rsv_data
#' @param model rsv_model
#' @param control rsv_control
#' @param subj_pre_list optional precomputed subject list (lambda_i, V_i)
#' @return numeric vector of length J (negative gradient)
grad_negloglik_beta <- function(beta,
                                eta,
                                data,
                                model,
                                control = make_rsv_control(),
                                subj_pre_list = NULL) {
  -grad_beta_fast(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    control       = control,
    subj_pre_list = subj_pre_list
  )
}


#' Gradient of the negative log-likelihood w.r.t. beta
#'
#' Computes the gradient of the **negative** total log-likelihood with respect
#' to the infection-age parameter vector \eqn{\beta}, holding \eqn{\eta} fixed.
#'
#' This is a lightweight wrapper around
#' \code{\link{grad_negloglik_beta}} intended for use with numerical
#' optimization routines such as \code{\link{optim}}. It assumes that all
#' subject-level precomputations (i.e., \code{lambda_i} and \code{V_i}) have
#' already been performed and are supplied via \code{subj_pre_list}.
#'
#' No additional validation or recomputation is performed beyond what is
#' handled internally by \code{grad_negloglik_beta()}.
#'
#' @param beta Numeric vector of length \eqn{J}. Infection-age parameter
#'   coefficients at which the gradient is evaluated.
#'
#' @param eta Numeric vector of length \eqn{K}. Visit-age parameter coefficients,
#'   treated as fixed in this gradient computation.
#'
#' @param data An \code{"rsv_data"} object containing the list of subjects.
#'
#' @param model An \code{"rsv_model"} object containing the global model
#'   components (basis matrices, circulation curve, etc.).
#'
#' @param rsv_control An \code{"rsv_control"} object governing numerical
#'   safeguards such as bound checking, clipping, and log-stabilization.
#'   This argument is named \code{rsv_control} to avoid clashes with the
#'   \code{control} argument used by \code{\link{optim}}.
#'
#' @param subj_pre_list A list of subject-level precomputations, typically
#'   produced by \code{\link{rsv_precompute_subject}}. Each element must contain
#'   at least:
#'   \itemize{
#'     \item \code{lambda_i} – subject-shifted circulation curve;
#'     \item \code{V_i}      – cumulative kernel matrix.
#'   }
#'   The list must have the same length and ordering as \code{data$subjects}.
#'
#' @return A numeric vector of length \eqn{J} giving the gradient of the
#'   **negative** total log-likelihood with respect to \eqn{\beta}.
#'
#' @details
#' This function is designed to be passed directly as the \code{gr} argument
#' to \code{\link{optim}} when optimizing over \eqn{\beta}:
#' \preformatted{
#' optim(par = beta_init,
#'       fn  = obj_negloglik_beta,
#'       gr  = obj_grad_negloglik_beta,
#'       ...)
#' }
#'
#' It is mathematically consistent with:
#' \itemize{
#'   \item \code{\link{obj_negloglik_beta}};
#'   \item \code{\link{grad_beta_fast}};
#'   \item \code{\link{rsv_negloglik}}.
#' }
#'
#' @seealso
#'   \code{\link{obj_negloglik_beta}},
#'   \code{\link{grad_beta_fast}},
#'   \code{\link{grad_negloglik_beta}},
#'   \code{\link{rsv_loglik}}
#'
#' @export
obj_grad_negloglik_beta <- function(beta,
                                    eta,
                                    data,
                                    model,
                                    rsv_control,
                                    subj_pre_list) {
  # --- basic checks (fail fast, cheap) ---
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  grad_negloglik_beta(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    control       = rsv_control,
    subj_pre_list = subj_pre_list
  )
}


#' Gradient of total log-likelihood with respect to eta (fast version)
#'
#' Computes the analytic gradient of the total (weighted) log-likelihood
#' \eqn{\ell(\beta,\eta)} with respect to the visit-age parameter vector
#' \eqn{\eta}, using a vectorized implementation for the no-visit (I2) case.
#'
#' @details
#' Visit case (I1):
#'   If subject i has a visit at age L_i,
#'   \deqn{
#'     \nabla_\eta \ell_i = s(L_i) / c(L_i;\eta),
#'   }
#'   where s(L_i) is the visit-age basis vector and
#'   c(L_i;\eta) = \eta^\top s(L_i).
#'
#' No-visit case (I2):
#'   If subject i has no visit in the first year of life,
#'   \deqn{
#'     \nabla_\eta \ell_i
#'     =
#'     -\frac{1}{\tau_i(\beta,\eta)}
#'     \sum_{m=1}^{365} s(m)\, Q_i(m;\beta),
#'   }
#'   where \eqn{\tau_i = 1 - \sum_m c(m;\eta)\, Q_i(m;\beta)}.
#'
#' The inner sum is evaluated efficiently as
#'   \code{S_day \%*\% Q},
#' where \code{S_day} is the K x 365 visit-age basis matrix.
#'
#' @param beta Numeric vector of length J. Infection-age basis coefficients
#'   (held fixed in this gradient).
#' @param eta Numeric vector of length K. Visit-age basis coefficients.
#' @param data An \code{"rsv_data"} object containing the subject list.
#' @param model An \code{"rsv_model"} object containing basis matrices.
#' @param control An \code{"rsv_control"} object controlling numerical checks.
#' @param subj_pre_list Optional list of subject-level precomputations.
#'
#' @return Numeric vector of length K giving \eqn{\nabla_\eta \ell(\beta,\eta)}.
#'
#' @export
grad_eta_fast <- function(beta,
                          eta,
                          data,
                          model,
                          control = make_rsv_control(),
                          subj_pre_list = NULL) {
  # --- basic checks ---
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.numeric(eta),  all(is.finite(eta)))
  stopifnot(!is.null(data$subjects), is.list(data$subjects))
  stopifnot(inherits(model, "rsv_model"))
  
  check_bounds <- isTRUE(control$check_bounds)
  warn_on_clip <- isTRUE(control$warn_on_clip)
  eps_log      <- if (!is.null(control$eps_log)) control$eps_log else 1e-12
  eps_tau      <- if (!is.null(control$eps_tau)) control$eps_tau else 1e-12
  
  subjects <- data$subjects
  n_subj   <- length(subjects)
  K        <- length(eta)
  
  # --- global precompute: w_vec, c_vec ---
  glob  <- rsv_precompute_global(beta = beta, eta = eta, model = model, control = control)
  w_vec <- glob$w_vec   # needed for Q in I2
  c_vec <- glob$c_vec   # stabilized visit weights
  
  # --- subject precompute if needed ---
  if (is.null(subj_pre_list)) {
    subj_pre_list <- lapply(subjects, rsv_precompute_subject,
                            model = model, control = control)
  } else {
    stopifnot(is.list(subj_pre_list), length(subj_pre_list) == n_subj)
  }
  
  grad  <- rep(0.0, K)
  S_day <- model$S_day  # K x 365
  
  for (i in seq_len(n_subj)) {
    subj   <- subjects[[i]]
    pre    <- subj_pre_list[[i]]
    weight <- subj$weight
    
    if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight)) {
      stop("subject weight must be finite scalar.")
    }
    
    lambda_i <- pre$lambda_i   # length 365
    V_i      <- pre$V_i        # J x 366
    
    L_i <- subj$visit_age
    
    if (!is.na(L_i)) {
      # -------------------------
      # I1 term: s(L) / c(L)
      # -------------------------
      L_i <- as.integer(L_i)
      
      s_L <- S_day[, L_i]
      c_L <- c_vec[L_i]
      
      if (!is.finite(c_L) || c_L <= 0) {
        if (check_bounds) {
          stop("c_vec[L_i] nonpositive/nonfinite in grad_eta_fast I1.")
        }
        if (warn_on_clip) {
          warning("c_vec[L_i] nonpositive/nonfinite in grad_eta_fast I1; repairing.")
        }
        c_L <- eps_log
      }
      
      grad <- grad + weight * (s_L / c_L)
      
    } else {
      # -------------------------
      # I2 term (vectorized)
      # -------------------------
      
      # Compute Fbar, pi, Q (same helpers as likelihood)
      Fbar <- Fbar_i(
        beta         = beta,
        V_i          = V_i,
        include_day0 = TRUE,
        check_bounds = check_bounds,
        warn_on_clip = warn_on_clip
      )
      
      pi <- pi_from_w(
        lambda_shift_i = lambda_i,
        w              = w_vec,
        check_bounds   = check_bounds,
        warn_on_clip   = warn_on_clip
      )
      
      Q <- Q_i(
        Fbar_i       = Fbar,
        pi_i         = pi,
        check_bounds = check_bounds,
        warn_on_clip = warn_on_clip
      )
      
      # g_i = sum_m c(m) Q(m)
      g_i <- sum(c_vec * Q)
      
      tau <- 1 - g_i
      tau <- .clip_to_range(
        tau, lower = eps_tau, upper = 1.0, name = "tau_i",
        check_bounds = check_bounds, warn_on_clip = warn_on_clip
      )
      
      # U_i(beta) = sum_m s(m) Q(m)
      # (Could be precomputed for eta-only optimization in a future version.)
      U <- as.vector(S_day %*% Q)
      
      grad <- grad + weight * (-(1 / tau) * U)
    }
  }
  
  names(grad) <- paste0("e", sprintf("%02d", seq_len(K)))
  
  grad
}


#' Gradient of negative log-likelihood w.r.t. eta
#'
#' Wrapper around grad_eta_fast() that returns the gradient of
#' the *negative* log-likelihood.
#'
#' @param beta numeric vector of length J (held fixed)
#' @param eta  numeric vector of length K
#' @param data rsv_data
#' @param model rsv_model
#' @param control rsv_control
#' @param subj_pre_list optional precomputed subject list
#'
#' @return numeric vector of length K (negative gradient)
#' @export
grad_negloglik_eta <- function(beta,
                               eta,
                               data,
                               model,
                               control = make_rsv_control(),
                               subj_pre_list = NULL) {
  -grad_eta_fast(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    control       = control,
    subj_pre_list = subj_pre_list
  )
}


#' Gradient of the negative log-likelihood w.r.t. eta
#'
#' Computes the gradient of the **negative** total log-likelihood
#' with respect to the visit-age parameter vector \eqn{\eta},
#' holding \eqn{\beta} fixed.
#'
#' @param eta Numeric vector of length K.
#' @param beta Numeric vector of length J (held fixed).
#' @param data rsv_data
#' @param model rsv_model
#' @param rsv_control rsv_control object
#' @param subj_pre_list list of subject-level precomputations
#'
#' @return numeric vector of length K
#' @export
obj_grad_negloglik_eta <- function(eta,
                                   beta,
                                   data,
                                   model,
                                   rsv_control,
                                   subj_pre_list) {
  stopifnot(is.numeric(eta), all(is.finite(eta)))
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  grad_negloglik_eta(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    control       = rsv_control,
    subj_pre_list = subj_pre_list
  )
}


#' Gradient of quadratic penalty: alpha * M theta
#'
#' @param theta Numeric parameter vector.
#' @param M Symmetric penalty matrix (typically t(D) %*% D).
#' @param alpha Nonnegative smoothing parameter.
#'
#' @return Numeric vector of same length as theta.
#' @export
penalty_grad <- function(theta, M, alpha) {
  # Fast exit
  if (alpha == 0) {
    return(rep(0.0, length(theta)))
  }
  
  # Minimal checks
  stopifnot(is.numeric(theta), is.numeric(alpha), alpha >= 0)
  stopifnot(is.matrix(M), ncol(M) == length(theta), nrow(M) == length(theta))
  
  # Compute alpha * M theta
  grad <- alpha * as.vector(M %*% theta)
  
  grad
}


#' Gradient of the penalized negative log-likelihood w.r.t. beta
#'
#' Computes the gradient of the **penalized negative log-likelihood**
#' with respect to the infection-age parameter vector \eqn{\beta},
#' holding \eqn{\eta} fixed.
#'
#' @details
#' The gradient corresponds to
#' \deqn{
#'   \nabla_\beta \tilde L(\beta)
#'   =
#'   -\nabla_\beta \ell(\beta,\eta)
#'   +
#'   \alpha_\beta M_\beta \beta,
#' }
#' where \eqn{M_\beta} is the second-difference penalty matrix stored
#' in \code{model$M_beta}.
#'
#' @param beta Numeric vector of length \eqn{J}.
#' @param eta Numeric vector of length \eqn{K} (held fixed).
#' @param alpha_beta Nonnegative scalar smoothing parameter.
#' @param data rsv_data
#' @param model rsv_model containing \code{M_beta}
#' @param rsv_control rsv_control object
#' @param subj_pre_list list of subject-level precomputations
#'
#' @return Numeric vector of length \eqn{J}.
#'
#' @seealso
#'   \code{\link{obj_pen_negloglik_beta}},
#'   \code{\link{grad_negloglik_beta}},
#'   \code{\link{penalty_grad}}
#'
#' @export
obj_grad_pen_negloglik_beta <- function(beta,
                                        eta,
                                        alpha_beta,
                                        data,
                                        model,
                                        rsv_control,
                                        subj_pre_list) {
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.numeric(alpha_beta), length(alpha_beta) == 1L, alpha_beta >= 0)
  stopifnot(!is.null(model$M_beta))
  
  # Gradient of negative log-likelihood
  g_nll <- obj_grad_negloglik_beta(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  # Gradient of penalty
  g_pen <- penalty_grad(
    theta = beta,
    M     = model$M_beta,
    alpha = alpha_beta
  )
  
  g_nll + g_pen
}


#' Gradient of the penalized negative log-likelihood w.r.t. eta
#'
#' Computes the gradient of the **penalized negative log-likelihood**
#' with respect to the visit-age parameter vector \eqn{\eta},
#' holding \eqn{\beta} fixed.
#'
#' @details
#' The gradient corresponds to
#' \deqn{
#'   \nabla_\eta \tilde L(\eta)
#'   =
#'   -\nabla_\eta \ell(\beta,\eta)
#'   +
#'   \alpha_\eta M_\eta \eta,
#' }
#' where:
#' \itemize{
#'   \item \eqn{\ell(\beta,\eta)} is the total weighted log-likelihood;
#'   \item \eqn{M_\eta} is the second-difference penalty matrix stored
#'         in \code{model$M_eta};
#'   \item \eqn{\alpha_\eta \ge 0} is the smoothing parameter.
#' }
#'
#' The likelihood gradient component is computed using
#' \code{\link{grad_negloglik_eta}}, ensuring full consistency with
#' the penalized objective.
#'
#' @param eta Numeric vector of length \eqn{K}.
#' @param beta Numeric vector of length \eqn{J} (held fixed).
#' @param alpha_eta Nonnegative scalar smoothing parameter.
#' @param data rsv_data
#' @param model rsv_model
#' @param rsv_control rsv_control object
#' @param subj_pre_list list of subject-level precomputations
#'
#' @return Numeric vector of length \eqn{K} giving the gradient of the
#'   penalized negative log-likelihood.
#'
#' @seealso
#'   \code{\link{obj_pen_negloglik_eta}},
#'   \code{\link{grad_negloglik_eta}},
#'   \code{\link{penalty_grad}}
#'
#' @export
obj_grad_pen_negloglik_eta <- function(eta,
                                       beta,
                                       alpha_eta,
                                       data,
                                       model,
                                       rsv_control,
                                       subj_pre_list) {
  stopifnot(is.numeric(eta), all(is.finite(eta)))
  stopifnot(is.numeric(alpha_eta), length(alpha_eta) == 1L, alpha_eta >= 0)
  stopifnot(!is.null(model$M_eta))
  
  # Gradient of negative log-likelihood
  g_nll <- obj_grad_negloglik_eta(
    eta           = eta,
    beta          = beta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  # Gradient of penalty
  g_pen <- penalty_grad(
    theta = eta,
    M     = model$M_eta,
    alpha = alpha_eta
  )
  
  g_nll + g_pen
}


#' Estimate eta via penalized maximum likelihood
#'
#' Optimizes the (possibly penalized) negative log-likelihood with respect
#' to the visit-age parameter vector eta, holding beta fixed.
#'
#' Always uses the analytic gradient and L-BFGS-B with nonnegativity
#' constraints (eta >= 0). The unpenalized case corresponds to
#' alpha_eta = 0.
#'
#' @param eta_init Numeric vector of length K. Initial value for eta.
#' @param beta Numeric vector of length J. Infection-age parameters (fixed).
#' @param alpha_eta Nonnegative scalar smoothing parameter.
#' @param data rsv_data object.
#' @param model rsv_model object (must contain M_eta).
#' @param rsv_control rsv_control object for numerical safeguards.
#' @param subj_pre_list List of subject-level precomputations.
#' @param maxit Maximum number of L-BFGS-B iterations.
#' @param trace Trace level passed to optim (0 = silent).
#' @param compute_nll Logical; if TRUE compute unpenalized NLL at solution.
#' @param fail_nonconvergence Logical; if TRUE stop() when convergence != 0.
#'
#' @return A list with components:
#'   \itemize{
#'     \item eta_hat: Estimated parameter vector.
#'     \item nll_pen: Penalized negative log-likelihood at solution.
#'     \item nll: Unpenalized negative log-likelihood (optional).
#'     \item convergence: Optim convergence code.
#'     \item message: Optim termination message.
#'     \item counts: Function/gradient evaluation counts.
#'   }
#'
#' @export
estimate_eta <- function(eta_init,
                         beta,
                         alpha_eta = 0,
                         data,
                         model,
                         rsv_control,
                         subj_pre_list,
                         maxit = 100,
                         trace = 0,
                         compute_nll = FALSE,
                         fail_nonconvergence = FALSE) {
  
  # --- basic checks ---
  stopifnot(is.numeric(eta_init), all(is.finite(eta_init)))
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.numeric(alpha_eta), length(alpha_eta) == 1L, alpha_eta >= 0)
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  K <- length(eta_init)
  
  # --- bounds: eta >= 0 ---
  lower <- rep(0, K)
  upper <- rep(Inf, K)
  
  # --- optimization ---
  fit <- optim(
    par    = eta_init,
    fn     = obj_pen_negloglik_eta,
    gr     = obj_grad_pen_negloglik_eta,
    method = "L-BFGS-B",
    lower  = lower,
    upper  = upper,
    beta   = beta,
    alpha_eta = alpha_eta,
    data   = data,
    model  = model,
    rsv_control = rsv_control,
    subj_pre_list = subj_pre_list,
    control = list(trace = trace, maxit = maxit)
  )
  
  # --- optional strict failure ---
  if (isTRUE(fail_nonconvergence) && fit$convergence != 0) {
    stop(
      sprintf(
        "estimate_eta(): optim did not converge (code=%d). Message: %s",
        fit$convergence,
        fit$message
      ),
      call. = FALSE
    )
  }
  
  # --- optional unpenalized NLL at solution ---
  value_nll <- NA_real_
  if (isTRUE(compute_nll)) {
    value_nll <- obj_negloglik_eta(
      eta           = fit$par,
      beta          = beta,
      data          = data,
      model         = model,
      rsv_control   = rsv_control,
      subj_pre_list = subj_pre_list
    )
  }
  
  # --- return structured result ---
  list(
    eta_hat     = fit$par,
    nll_pen     = fit$value,
    nll         = value_nll,
    convergence = fit$convergence,
    message     = fit$message,
    counts      = fit$counts
  )
}


#' Estimate beta via penalized maximum likelihood
#'
#' Optimizes the (possibly penalized) negative log-likelihood with respect
#' to the infection-age parameter vector beta, holding eta fixed.
#'
#' Always uses the analytic gradient and L-BFGS-B with nonnegativity
#' constraints (beta >= 0). The unpenalized case corresponds to
#' alpha_beta = 0.
#'
#' @param beta_init Numeric vector of length J. Initial value for beta.
#' @param eta Numeric vector of length K. Visit-age parameters (fixed).
#' @param alpha_beta Nonnegative scalar smoothing parameter.
#' @param data rsv_data object.
#' @param model rsv_model object (must contain M_beta).
#' @param rsv_control rsv_control object for numerical safeguards.
#' @param subj_pre_list List of subject-level precomputations.
#' @param maxit Maximum number of L-BFGS-B iterations.
#' @param trace Trace level passed to optim (0 = silent).
#' @param compute_nll Logical; if TRUE compute unpenalized NLL at solution.
#' @param fail_nonconvergence Logical; if TRUE stop() when convergence != 0.
#'
#' @return A list with components:
#'   \itemize{
#'     \item beta_hat: Estimated parameter vector.
#'     \item nll_pen: Penalized negative log-likelihood at solution.
#'     \item nll: Unpenalized negative log-likelihood (optional).
#'     \item convergence: Optim convergence code.
#'     \item message: Optim termination message.
#'     \item counts: Function/gradient evaluation counts.
#'   }
#'
#' @export
estimate_beta <- function(beta_init,
                          eta,
                          alpha_beta = 0,
                          data,
                          model,
                          rsv_control,
                          subj_pre_list,
                          maxit = 200,
                          trace = 0,
                          compute_nll = FALSE,
                          fail_nonconvergence = FALSE) {
  
  # --- basic checks ---
  stopifnot(is.numeric(beta_init), all(is.finite(beta_init)))
  stopifnot(is.numeric(eta), all(is.finite(eta)))
  stopifnot(is.numeric(alpha_beta), length(alpha_beta) == 1L, alpha_beta >= 0)
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  J <- length(beta_init)
  
  # --- bounds: beta >= 0 ---
  lower <- rep(0, J)
  upper <- rep(Inf, J)
  
  # --- optimization ---
  fit <- optim(
    par    = beta_init,
    fn     = obj_pen_negloglik_beta,
    gr     = obj_grad_pen_negloglik_beta,
    method = "L-BFGS-B",
    lower  = lower,
    upper  = upper,
    eta    = eta,
    alpha_beta = alpha_beta,
    data   = data,
    model  = model,
    rsv_control = rsv_control,
    subj_pre_list = subj_pre_list,
    control = list(trace = trace, maxit = maxit)
  )
  
  # --- optional strict failure ---
  if (isTRUE(fail_nonconvergence) && fit$convergence != 0) {
    stop(
      sprintf(
        "estimate_beta(): optim did not converge (code=%d). Message: %s",
        fit$convergence,
        fit$message
      ),
      call. = FALSE
    )
  }
  
  # --- optional unpenalized NLL at solution ---
  value_nll <- NA_real_
  if (isTRUE(compute_nll)) {
    value_nll <- obj_negloglik_beta(
      beta          = fit$par,
      eta           = eta,
      data          = data,
      model         = model,
      rsv_control   = rsv_control,
      subj_pre_list = subj_pre_list
    )
  }
  
  # --- return structured result ---
  list(
    beta_hat    = fit$par,
    nll_pen     = fit$value,
    nll         = value_nll,
    convergence = fit$convergence,
    message     = fit$message,
    counts      = fit$counts
  )
}


#' Alternating estimation of beta and eta via block coordinate descent
#'
#' Fits the joint penalized model for the infection-age spline
#' coefficients \eqn{\beta} and the covariate spline coefficients
#' \eqn{\eta} using alternating (block coordinate descent) optimization.
#'
#' At each outer iteration:
#' \enumerate{
#'   \item Update \eqn{\beta} by minimizing the penalized objective
#'         with \eqn{\eta} held fixed.
#'   \item Update \eqn{\eta} by minimizing the penalized objective
#'         with \eqn{\beta} held fixed.
#' }
#'
#' The joint penalized objective minimized is
#' \deqn{
#'   F(\beta, \eta)
#'   =
#'   \mathrm{NLL}(\beta, \eta)
#'   +
#'   \frac{\alpha}{2}
#'   \left(
#'     \beta^\top M_\beta \beta
#'     +
#'     \eta^\top M_\eta \eta
#'   \right),
#' }
#' where \eqn{\mathrm{NLL}(\beta, \eta)} is the unpenalized negative
#' log-likelihood returned by \code{\link{rsv_negloglik}}.
#'
#' Convergence may be assessed using relative parameter change
#' (L2 norm), relative objective change, or both.
#'
#' @param beta_init Numeric vector giving the initial values for
#'   the spline coefficients \eqn{\beta}.
#'
#' @param eta_init Numeric vector giving the initial values for
#'   the spline coefficients \eqn{\eta}.
#'
#' @param alpha Non-negative scalar smoothing parameter controlling
#'   the strength of penalization for both \eqn{\beta} and \eqn{\eta}.
#'
#' @param data Data object passed to \code{\link{rsv_negloglik}} and
#'   the blockwise estimation routines.
#'
#' @param model Model object containing spline basis matrices and
#'   penalty matrices \code{M_beta} and \code{M_eta}.
#'
#' @param rsv_control Control object passed to the likelihood and
#'   estimation routines.
#'
#' @param subj_pre_list Precomputed subject-level quantities used to
#'   accelerate likelihood evaluation.
#'
#' @param max_outer Maximum number of outer (alternating) iterations.
#'
#' @param convergence Character string specifying the convergence
#'   criterion. One of:
#'   \describe{
#'     \item{\code{"objective"}}{Stop when the relative change in the
#'       penalized objective is below \code{tol_obj}.}
#'     \item{\code{"l2"}}{Stop when the relative L2 change in both
#'       \eqn{\beta} and \eqn{\eta} is below \code{tol_l2}.}
#'     \item{\code{"both"}}{Require both objective and L2 criteria.}
#'   }
#'
#' @param tol_obj Relative tolerance for objective-based convergence.
#'
#' @param tol_l2 Relative tolerance for L2-based convergence.
#'
#' @param verbose Logical; if \code{TRUE}, prints iteration-level
#'   diagnostics.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{beta_hat}}{Estimated spline coefficients \eqn{\beta}.}
#'     \item{\code{eta_hat}}{Estimated spline coefficients \eqn{\eta}.}
#'     \item{\code{converged}}{Logical indicating whether convergence
#'       criteria were satisfied.}
#'     \item{\code{n_outer}}{Number of outer iterations performed.}
#'     \item{\code{inner_conv_beta}}{Logical vector indicating whether
#'       each beta-update optimization converged.}
#'     \item{\code{inner_conv_eta}}{Logical vector indicating whether
#'       each eta-update optimization converged.}
#'     \item{\code{inner_failed}}{Logical; \code{TRUE} if any inner
#'       optimization reported non-convergence.}
#'     \item{\code{rel_change_beta}}{Final relative L2 change in
#'       \eqn{\beta}.}
#'     \item{\code{rel_change_eta}}{Final relative L2 change in
#'       \eqn{\eta}.}
#'     \item{\code{rel_change_obj}}{Final relative change in the
#'       penalized objective.}
#'     \item{\code{final_objective}}{Value of the penalized objective
#'       at the final parameter estimates.}
#'   }
#'
#' @details
#' The algorithm implements block coordinate descent and is
#' typically monotone in the penalized objective, though small
#' numerical increases may occur due to optimizer tolerances.
#'
#' If a non-finite objective value is encountered, the routine
#' terminates early and returns the current parameter values with
#' \code{converged = FALSE}.
#'
#' @seealso \code{\link{estimate_beta}},
#'   \code{\link{estimate_eta}},
#'   \code{\link{full_penalized_objective}}
#'
#' @export
estimate_model_alternating <- function(
    beta_init,
    eta_init,
    alpha,
    data,
    model,
    rsv_control,
    subj_pre_list,
    max_outer = 25L,
    convergence = c("objective", "l2", "both"),
    tol_obj = 1e-5,
    tol_l2  = 1e-6,
    verbose = FALSE
) {
  convergence <- match.arg(convergence)
  
  stopifnot(is.numeric(beta_init))
  stopifnot(is.numeric(eta_init))
  stopifnot(is.list(subj_pre_list))
  
  eps <- 1e-12
  
  beta_curr <- beta_init
  eta_curr  <- eta_init
  
  # Compute initial objective
  obj_old <- full_penalized_objective(
    beta_curr, eta_curr, alpha,
    data, model, rsv_control, subj_pre_list
  )
  
  if (!is.finite(obj_old)) {
    return(list(
      beta_hat = beta_curr,
      eta_hat  = eta_curr,
      converged = FALSE,
      n_outer = 0L,
      inner_conv_beta = logical(0),
      inner_conv_eta  = logical(0),
      inner_failed = TRUE,
      rel_change_beta = NA_real_,
      rel_change_eta  = NA_real_,
      rel_change_obj  = NA_real_,
      final_objective = obj_old
    ))
  }
  
  inner_conv_beta <- logical(max_outer)
  inner_conv_eta  <- logical(max_outer)
  
  rel_change_beta <- NA_real_
  rel_change_eta  <- NA_real_
  rel_change_obj  <- NA_real_
  
  converged <- FALSE
  
  for (iter in seq_len(max_outer)) {
    
    if (verbose) {
      message(sprintf("Outer iteration %d", iter))
    }
    
    ## ---- β-step ----
    fit_beta <- estimate_beta(
      beta_init     = beta_curr,
      eta           = eta_curr,
      data          = data,
      model         = model,
      rsv_control   = rsv_control,
      subj_pre_list = subj_pre_list,
      alpha_beta    = alpha
    )
    
    beta_new <- fit_beta$beta_hat
    inner_conv_beta[iter] <- (fit_beta$convergence == 0)
    
    ## ---- η-step ----
    fit_eta <- estimate_eta(
      eta_init      = eta_curr,
      beta          = beta_new,
      data          = data,
      model         = model,
      rsv_control   = rsv_control,
      subj_pre_list = subj_pre_list,
      alpha_eta     = alpha
    )
    
    eta_new <- fit_eta$eta_hat
    inner_conv_eta[iter] <- (fit_eta$convergence == 0)
    
    ## ---- Compute full objective ----
    obj_new <- full_penalized_objective(
      beta_new, eta_new, alpha,
      data, model, rsv_control, subj_pre_list
    )
    
    if (!is.finite(obj_new)) {
      if (verbose) {
        message("Non-finite objective encountered. Stopping.")
      }
      n_outer <- iter
      break
    }
    
    ## ---- Relative L2 changes ----
    denom_beta <- sqrt(sum(beta_curr^2))
    if (denom_beta < eps) denom_beta <- 1
    
    denom_eta <- sqrt(sum(eta_curr^2))
    if (denom_eta < eps) denom_eta <- 1
    
    rel_change_beta <- sqrt(sum((beta_new - beta_curr)^2)) / denom_beta
    rel_change_eta  <- sqrt(sum((eta_new  - eta_curr)^2))  / denom_eta
    
    ## ---- Relative objective change ----
    denom_obj <- abs(obj_old)
    if (denom_obj < eps) denom_obj <- 1
    
    rel_change_obj <- abs(obj_new - obj_old) / denom_obj
    
    if (verbose) {
      message(sprintf(
        "  rel_obj = %.3e | rel_beta = %.3e | rel_eta = %.3e",
        rel_change_obj, rel_change_beta, rel_change_eta
      ))
    }
    
    ## ---- Update current values ----
    beta_curr <- beta_new
    eta_curr  <- eta_new
    obj_old   <- obj_new
    
    ## ---- Check convergence ----
    stop_l2  <- (rel_change_beta < tol_l2 && rel_change_eta < tol_l2)
    stop_obj <- (rel_change_obj  < tol_obj)
    
    should_stop <- switch(
      convergence,
      l2        = stop_l2,
      objective = stop_obj,
      both      = (stop_l2 && stop_obj)
    )
    
    if (should_stop) {
      converged <- TRUE
      n_outer <- iter
      break
    }
    
    if (iter == max_outer) {
      n_outer <- max_outer
    }
  }
  
  ## Trim diagnostic vectors
  inner_conv_beta <- inner_conv_beta[seq_len(n_outer)]
  inner_conv_eta  <- inner_conv_eta[seq_len(n_outer)]
  
  inner_failed <- any(!inner_conv_beta) || any(!inner_conv_eta)
  
  return(list(
    beta_hat         = beta_curr,
    eta_hat          = eta_curr,
    converged        = converged,
    n_outer          = n_outer,
    inner_conv_beta  = inner_conv_beta,
    inner_conv_eta   = inner_conv_eta,
    inner_failed     = inner_failed,
    rel_change_beta  = rel_change_beta,
    rel_change_eta   = rel_change_eta,
    rel_change_obj   = rel_change_obj,
    final_objective  = obj_old
  ))
}


#' Create perturbed initial coefficient values
#'
#' Generates reproducible starting values for alternating optimization by
#' applying independent multiplicative log-normal perturbations to two
#' nonnegative coefficient vectors. The perturbed vectors are rescaled to
#' preserve approximately the same mean coefficient level as the originals.
#'
#' This helper is intended for controlled simulation examples where the
#' optimization should begin near known reference or true parameter values.
#'
#' @param beta A nonempty numeric vector of nonnegative coefficients for the
#'   infection-age curve.
#' @param eta A nonempty numeric vector of nonnegative coefficients for the
#'   visit-age curve.
#' @param noise_sd A single nonnegative numeric value giving the standard
#'   deviation of the normal noise applied on the log scale. Defaults to
#'   `0.05`.
#' @param seed An optional numeric value used to initialize the random-number
#'   generator. Defaults to `NULL`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{beta_init}{The perturbed initial coefficients corresponding to
#'     `beta`.}
#'   \item{eta_init}{The perturbed initial coefficients corresponding to
#'     `eta`.}
#' }
#'
#' @details
#' For each coefficient vector \eqn{\theta}, the function constructs initial
#' values proportional to
#'
#' \deqn{\theta_j^{init} = \theta_j \exp(\epsilon_j),}
#'
#' where \eqn{\epsilon_j} are independent normal random variables with mean
#' zero and standard deviation `noise_sd`. This multiplicative construction
#' preserves nonnegativity and produces perturbations relative to the scale of
#' each coefficient.
#'
#' Coefficients that are zero or numerically close to zero are bounded below
#' before the perturbed vector is rescaled to match the original mean.
#'
#' @examples
#' initial_values <- make_perturbed_initial_values(
#'   beta = c(0.1, 0.2, 0.3),
#'   eta = c(0.3, 0.2, 0.1),
#'   noise_sd = 0.05,
#'   seed = 123
#' )
#'
#' initial_values$beta_init
#' initial_values$eta_init
#'
#' @keywords internal
make_perturbed_initial_values <- function(
    beta,
    eta,
    noise_sd = 0.05,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  perturb <- function(coef) {
    coef_init <- coef * exp(
      rnorm(
        n    = length(coef),
        mean = 0,
        sd   = noise_sd
      )
    )
    
    coef_init <- pmax(coef_init, 1e-8)
    
    coef_init * mean(coef) / mean(coef_init)
  }
  
  list(
    beta_init = perturb(beta),
    eta_init  = perturb(eta)
  )
}
