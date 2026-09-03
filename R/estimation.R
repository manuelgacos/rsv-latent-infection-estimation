# Estimation
#
# Implements analytic gradients, penalized objective gradients, blockwise
# optimization, alternating estimation, and initialization for the RSV model.


# Likelihood gradients -----------------------------------------------------


#' Compute the log-likelihood gradient with respect to beta
#'
#' Computes the analytic gradient of the total weighted log-likelihood with
#' respect to the infection-age spline coefficients, holding the
#' healthcare-visit coefficients fixed. The no-visit contribution is
#' evaluated using matrix-vector products for efficient repeated estimation.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients, held fixed in this gradient.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases and
#'   seasonal RSV circulation curve.
#' @param control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list Optional list of subject-level precomputations, with
#'   one element per subject. If `NULL`, the precomputations are constructed
#'   internally.
#'
#' @return Numeric vector of length `J` containing the gradient of the total
#'   weighted log-likelihood with respect to \eqn{\beta}.
#'
#' @details
#' The gradient is the weighted sum of subject-level contributions. For a
#' subject with a healthcare visit at age \eqn{L_i},
#' \deqn{
#'   \nabla_\beta \ell_i
#'   =
#'   \frac{b(L_i)}{w(L_i;\beta)} - v_i(L_i).
#' }
#'
#' For a subject with no healthcare visit during the first year of life,
#' \deqn{
#'   \nabla_\beta \ell_i
#'   =
#'   -\frac{1}{\tau_i(\beta,\eta)}
#'   \sum_{m=1}^{365}
#'   c(m;\eta)
#'   \left[
#'     \bar F_i(m;\beta)\lambda_i(m)b(m)
#'     -
#'     Q_i(m;\beta)v_i(m-1)
#'   \right].
#' }
#'
#' The no-visit sum is evaluated using matrix-vector products involving
#' `model$B_day` and the first 365 columns of `V_i`. The daily curves use the
#' numerical safeguards applied by `rsv_precompute_global()`, and
#' \eqn{\tau_i} is bounded below by `control$eps_tau`. Numerical clipping is
#' treated as fixed when forming the analytic gradient.
grad_beta_fast <- function(beta,
                           eta,
                           data,
                           model,
                           control = make_rsv_control(),
                           subj_pre_list = NULL) {
  
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
  
  # Precompute coefficient-dependent curves once for all subjects.
  glob  <- rsv_precompute_global(beta = beta, eta = eta, model = model, control = control)
  w_vec <- glob$w_vec
  c_vec <- glob$c_vec
  
  # Reuse supplied subject precomputations or construct them once.
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
    
    lambda_i <- pre$lambda_i        
    V_i      <- pre$V_i             
    
    L_i <- subj$visit_age
    
    if (!is.na(L_i)) {
      # I1: visit contribution b(L) / w(L) - v(L).
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
      # I2: evaluate the no-visit contribution with matrix-vector products.
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
      
      # Bound tau away from zero before division in the no-visit gradient.
      tau <- .clip_to_range(
        tau, lower = eps_tau, upper = 1.0, name = "tau_i",
        check_bounds = check_bounds, warn_on_clip = warn_on_clip
      )
      
      # a[m] = c(m) * Fbar(m) * lambda_i(m), with
      # Fbar(m) = Fbar[m + 1].
      a <- c_vec * (Fbar[2:366] * lambda_i)
      
      # b[m] = c(m) * Q(m).
      b <- c_vec * Q                        
      
      # Combine the two J-dimensional terms in the no-visit gradient sum.
      inner <- as.vector(B_day %*% a) - as.vector(V_i[, 1:365, drop = FALSE] %*% b)
      
      grad <- grad + weight * (-(1 / tau)) * inner
    }
  }
  
  grad
}


#' Compute the negative log-likelihood gradient with respect to beta
#'
#' Returns the gradient of the negative total weighted log-likelihood with
#' respect to the infection-age spline coefficients. This is the
#' sign-reversed gradient returned by `grad_beta_fast()`.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients, held fixed in this gradient.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases and
#'   seasonal RSV circulation curve.
#' @param control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list Optional list of subject-level precomputations, with
#'   one element per subject. If `NULL`, the precomputations are constructed
#'   internally by `grad_beta_fast()`.
#'
#' @return Numeric vector of length `J` containing the gradient of the
#'   negative total weighted log-likelihood with respect to \eqn{\beta}.
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


#' Compute the beta gradient for optimization
#'
#' Evaluates the gradient of the negative total weighted log-likelihood with
#' respect to the infection-age spline coefficients. This optimizer-facing
#' wrapper delegates the calculation to `grad_negloglik_beta()`.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients, held fixed in this gradient.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases and
#'   seasonal RSV circulation curve.
#' @param rsv_control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization. The argument name avoids conflict with the
#'   `control` argument used by `optim()`.
#' @param subj_pre_list List of subject-level precomputations, with one element
#'   per subject.
#'
#' @return Numeric vector of length `J` containing the gradient of the
#'   negative total weighted log-likelihood with respect to \eqn{\beta}.
obj_grad_negloglik_beta <- function(beta,
                                    eta,
                                    data,
                                    model,
                                    rsv_control,
                                    subj_pre_list) {
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


#' Compute the log-likelihood gradient with respect to eta
#'
#' Computes the analytic gradient of the total weighted log-likelihood with
#' respect to the healthcare-visit spline coefficients, holding the
#' infection-age coefficients fixed. The no-visit contribution is evaluated
#' using a matrix-vector product for efficient repeated estimation.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients, held fixed in this gradient.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases and
#'   seasonal RSV circulation curve.
#' @param control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list Optional list of subject-level precomputations, with
#'   one element per subject. If `NULL`, the precomputations are constructed
#'   internally.
#'
#' @return Numeric vector of length `K` containing the gradient of the total
#'   weighted log-likelihood with respect to \eqn{\eta}.
#'
#' @details
#' The gradient is the weighted sum of subject-level contributions. For a
#' subject with a healthcare visit at age \eqn{L_i},
#' \deqn{
#'   \nabla_\eta \ell_i
#'   =
#'   \frac{s(L_i)}{c(L_i;\eta)}.
#' }
#'
#' For a subject with no healthcare visit during the first year of life,
#' \deqn{
#'   \nabla_\eta \ell_i
#'   =
#'   -\frac{1}{\tau_i(\beta,\eta)}
#'   \sum_{m=1}^{365}
#'   s(m)Q_i(m;\beta).
#' }
#'
#' The no-visit sum is evaluated as `model$S_day %*% Q`. The daily
#' healthcare-visit curve uses the numerical safeguards applied by
#' `rsv_precompute_global()`, and \eqn{\tau_i} is bounded below by
#' `control$eps_tau`.
grad_eta_fast <- function(beta,
                          eta,
                          data,
                          model,
                          control = make_rsv_control(),
                          subj_pre_list = NULL) {
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
  
  # Precompute coefficient-dependent curves once for all subjects.
  glob  <- rsv_precompute_global(beta = beta, eta = eta, model = model, control = control)
  w_vec <- glob$w_vec   
  c_vec <- glob$c_vec   
  
  # Reuse supplied subject precomputations or construct them once.
  if (is.null(subj_pre_list)) {
    subj_pre_list <- lapply(subjects, rsv_precompute_subject,
                            model = model, control = control)
  } else {
    stopifnot(is.list(subj_pre_list), length(subj_pre_list) == n_subj)
  }
  
  grad  <- rep(0.0, K)
  # Healthcare-visit basis with K rows and 365 age-day columns.
  S_day <- model$S_day
  
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
      # I1: visit contribution s(L) / c(L).
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
      # I2: evaluate the no-visit contribution with a matrix-vector product.
      
      # Reuse the likelihood quantities needed for the no-visit gradient.
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
      
      # g_i = sum_m c(m) Q(m).
      g_i <- sum(c_vec * Q)
      
      tau <- 1 - g_i
      # Bound tau away from zero before division in the no-visit gradient.
      tau <- .clip_to_range(
        tau, lower = eps_tau, upper = 1.0, name = "tau_i",
        check_bounds = check_bounds, warn_on_clip = warn_on_clip
      )
      
      # U_i(beta) = sum_m s(m) * Q(m).
      U <- as.vector(S_day %*% Q)
      
      grad <- grad + weight * (-(1 / tau) * U)
    }
  }
  
  names(grad) <- paste0("e", sprintf("%02d", seq_len(K)))
  
  grad
}


#' Compute the negative log-likelihood gradient with respect to eta
#'
#' Returns the gradient of the negative total weighted log-likelihood with
#' respect to the healthcare-visit spline coefficients. This is the
#' sign-reversed gradient returned by `grad_eta_fast()`.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients, held fixed in this gradient.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases and
#'   seasonal RSV circulation curve.
#' @param control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list Optional list of subject-level precomputations, with
#'   one element per subject. If `NULL`, the precomputations are constructed
#'   internally by `grad_eta_fast()`.
#'
#' @return Numeric vector of length `K` containing the gradient of the
#'   negative total weighted log-likelihood with respect to \eqn{\eta}.
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


#' Compute the eta gradient for optimization
#'
#' Evaluates the gradient of the negative total weighted log-likelihood with
#' respect to the healthcare-visit spline coefficients. This optimizer-facing
#' wrapper delegates the calculation to `grad_negloglik_eta()`.
#'
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients, held fixed in this gradient.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases and
#'   seasonal RSV circulation curve.
#' @param rsv_control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization. The argument name avoids conflict with the
#'   `control` argument used by `optim()`.
#' @param subj_pre_list List of subject-level precomputations, with one element
#'   per subject.
#'
#' @return Numeric vector of length `K` containing the gradient of the
#'   negative total weighted log-likelihood with respect to \eqn{\eta}.
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


# Penalty gradients --------------------------------------------------------


#' Compute the quadratic penalty gradient
#'
#' Computes the gradient of the quadratic smoothing penalty with respect to a
#' coefficient vector.
#'
#' @param theta Numeric vector of length `p` containing the coefficients at
#'   which the penalty gradient is evaluated.
#' @param M Numeric `p x p` penalty matrix, typically constructed as
#'   \eqn{D^\top D} from a second-difference matrix.
#' @param alpha Nonnegative numeric scalar giving the smoothing parameter.
#'
#' @return Numeric vector of length `p` containing
#'   \eqn{\alpha M\theta}. If `alpha` is zero, returns a zero vector.
penalty_grad <- function(theta, M, alpha) {
  # Return the zero gradient when the penalty is inactive.
  if (alpha == 0) {
    return(rep(0.0, length(theta)))
  }
  
  stopifnot(is.numeric(theta), is.numeric(alpha), alpha >= 0)
  stopifnot(is.matrix(M), ncol(M) == length(theta), nrow(M) == length(theta))
  
  grad <- alpha * as.vector(M %*% theta)
  
  grad
}


#' Compute the penalized gradient for beta optimization
#'
#' Computes the gradient of the penalized negative log-likelihood with respect
#' to the infection-age spline coefficients, holding the healthcare-visit
#' coefficients fixed.
#'
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients, held fixed during optimization.
#' @param alpha_beta Nonnegative numeric scalar controlling the smoothing
#'   penalty for \eqn{\beta}.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the penalty matrix `M_beta`
#'   and the components used in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one element
#'   per subject.
#'
#' @return Numeric vector of length `J` containing the gradient of the
#'   penalized negative log-likelihood with respect to \eqn{\beta}.
#'
#' @details
#' The gradient is
#' \deqn{
#'   -\nabla_\beta \ell(\beta,\eta)
#'   +
#'   \alpha_\beta M_\beta \beta,
#' }
#' where \eqn{M_\beta} is the quadratic smoothing-penalty matrix stored in
#' `model$M_beta`.
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
  
  g_nll <- obj_grad_negloglik_beta(
    beta          = beta,
    eta           = eta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  g_pen <- penalty_grad(
    theta = beta,
    M     = model$M_beta,
    alpha = alpha_beta
  )
  
  g_nll + g_pen
}


#' Compute the penalized gradient for eta optimization
#'
#' Computes the gradient of the penalized negative log-likelihood with respect
#' to the healthcare-visit spline coefficients, holding the infection-age
#' coefficients fixed.
#'
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients, held fixed during optimization.
#' @param alpha_eta Nonnegative numeric scalar controlling the smoothing
#'   penalty for \eqn{\eta}.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the penalty matrix `M_eta`
#'   and the components used in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one element
#'   per subject.
#'
#' @return Numeric vector of length `K` containing the gradient of the
#'   penalized negative log-likelihood with respect to \eqn{\eta}.
#'
#' @details
#' The gradient is
#' \deqn{
#'   -\nabla_\eta \ell(\beta,\eta)
#'   +
#'   \alpha_\eta M_\eta \eta,
#' }
#' where \eqn{M_\eta} is the quadratic smoothing-penalty matrix stored in
#' `model$M_eta`.
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
  
  g_nll <- obj_grad_negloglik_eta(
    eta           = eta,
    beta          = beta,
    data          = data,
    model         = model,
    rsv_control   = rsv_control,
    subj_pre_list = subj_pre_list
  )
  
  g_pen <- penalty_grad(
    theta = eta,
    M     = model$M_eta,
    alpha = alpha_eta
  )
  
  g_nll + g_pen
}


# Blockwise estimation -----------------------------------------------------


#' Estimate the healthcare-visit coefficients
#'
#' Estimates the healthcare-visit spline coefficients by minimizing the
#' penalized negative log-likelihood while holding the infection-age
#' coefficients fixed. Optimization uses the analytic gradient with L-BFGS-B
#' and constrains the healthcare-visit coefficients to be nonnegative.
#'
#' @param eta_init Numeric vector of length `K` containing the initial
#'   healthcare-visit spline coefficients. Values should be nonnegative.
#' @param beta Numeric vector of length `J` containing the infection-age
#'   spline coefficients, held fixed during optimization.
#' @param alpha_eta Nonnegative numeric scalar controlling the smoothing
#'   penalty for \eqn{\eta}.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the healthcare-visit basis,
#'   penalty matrix `M_eta`, and components used in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one element
#'   per subject.
#' @param maxit Positive integer scalar giving the maximum number of L-BFGS-B
#'   iterations.
#' @param trace Nonnegative integer scalar giving the trace level passed to
#'   `optim()`. A value of zero suppresses optimizer tracing.
#' @param compute_nll Logical scalar indicating whether to evaluate the
#'   unpenalized negative log-likelihood at the estimated coefficients.
#' @param fail_nonconvergence Logical scalar indicating whether to stop when
#'   `optim()` returns a nonzero convergence code.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{eta_hat}}{Numeric vector of length `K` containing the
#'     estimated healthcare-visit spline coefficients.}
#'   \item{\code{nll_pen}}{Numeric scalar giving the penalized negative
#'     log-likelihood at `eta_hat`.}
#'   \item{\code{nll}}{Numeric scalar giving the unpenalized negative
#'     log-likelihood at `eta_hat` when `compute_nll = TRUE`; otherwise
#'     `NA_real_`.}
#'   \item{\code{convergence}}{Integer scalar containing the `optim()`
#'     convergence code.}
#'   \item{\code{message}}{Termination message returned by `optim()`.}
#'   \item{\code{counts}}{Named vector containing the `optim()` function and
#'     gradient evaluation counts.}
#' }
#'
#' @details
#' With \eqn{\beta} fixed, the function minimizes
#' \deqn{
#'   \mathrm{NLL}(\beta,\eta)
#'   +
#'   \frac{\alpha_\eta}{2}
#'   \eta^\top M_\eta \eta
#' }
#' subject to \eqn{\eta \ge 0}. Setting `alpha_eta = 0` gives the unpenalized
#' blockwise estimate.
#'
#' If `fail_nonconvergence = FALSE`, a nonzero optimizer convergence code does
#' not stop execution and is instead returned with the optimizer message.
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
  
  stopifnot(is.numeric(eta_init), all(is.finite(eta_init)))
  stopifnot(is.numeric(beta), all(is.finite(beta)))
  stopifnot(is.numeric(alpha_eta), length(alpha_eta) == 1L, alpha_eta >= 0)
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  K <- length(eta_init)
  
  # Enforce nonnegativity of the healthcare-visit coefficients.
  lower <- rep(0, K)
  upper <- rep(Inf, K)
  
  # Optimize eta with beta fixed using the analytic penalized gradient.
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
  
  # Optionally treat optimizer nonconvergence as an error.
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
  
  # Optionally evaluate the unpenalized NLL for diagnostics.
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
  
  list(
    eta_hat     = fit$par,
    nll_pen     = fit$value,
    nll         = value_nll,
    convergence = fit$convergence,
    message     = fit$message,
    counts      = fit$counts
  )
}


#' Estimate the infection-age coefficients
#'
#' Estimates the infection-age spline coefficients by minimizing the
#' penalized negative log-likelihood while holding the healthcare-visit
#' coefficients fixed. Optimization uses the analytic gradient with L-BFGS-B
#' and constrains the infection-age coefficients to be nonnegative.
#'
#' @param beta_init Numeric vector of length `J` containing the initial
#'   infection-age spline coefficients. Values should be nonnegative.
#' @param eta Numeric vector of length `K` containing the healthcare-visit
#'   spline coefficients, held fixed during optimization.
#' @param alpha_beta Nonnegative numeric scalar controlling the smoothing
#'   penalty for \eqn{\beta}.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the infection-age basis,
#'   penalty matrix `M_beta`, and components used in likelihood evaluation.
#' @param rsv_control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one element
#'   per subject.
#' @param maxit Positive integer scalar giving the maximum number of L-BFGS-B
#'   iterations.
#' @param trace Nonnegative integer scalar giving the trace level passed to
#'   `optim()`. A value of zero suppresses optimizer tracing.
#' @param compute_nll Logical scalar indicating whether to evaluate the
#'   unpenalized negative log-likelihood at the estimated coefficients.
#' @param fail_nonconvergence Logical scalar indicating whether to stop when
#'   `optim()` returns a nonzero convergence code.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{beta_hat}}{Numeric vector of length `J` containing the
#'     estimated infection-age spline coefficients.}
#'   \item{\code{nll_pen}}{Numeric scalar giving the penalized negative
#'     log-likelihood at `beta_hat`.}
#'   \item{\code{nll}}{Numeric scalar giving the unpenalized negative
#'     log-likelihood at `beta_hat` when `compute_nll = TRUE`; otherwise
#'     `NA_real_`.}
#'   \item{\code{convergence}}{Integer scalar containing the `optim()`
#'     convergence code.}
#'   \item{\code{message}}{Termination message returned by `optim()`.}
#'   \item{\code{counts}}{Named vector containing the `optim()` function and
#'     gradient evaluation counts.}
#' }
#'
#' @details
#' With \eqn{\eta} fixed, the function minimizes
#' \deqn{
#'   \mathrm{NLL}(\beta,\eta)
#'   +
#'   \frac{\alpha_\beta}{2}
#'   \beta^\top M_\beta \beta
#' }
#' subject to \eqn{\beta \ge 0}. Setting `alpha_beta = 0` gives the
#' unpenalized blockwise estimate.
#'
#' If `fail_nonconvergence = FALSE`, a nonzero optimizer convergence code does
#' not stop execution and is instead returned with the optimizer message.
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
  
  stopifnot(is.numeric(beta_init), all(is.finite(beta_init)))
  stopifnot(is.numeric(eta), all(is.finite(eta)))
  stopifnot(is.numeric(alpha_beta), length(alpha_beta) == 1L, alpha_beta >= 0)
  stopifnot(is.list(subj_pre_list))
  stopifnot(length(subj_pre_list) == length(data$subjects))
  
  J <- length(beta_init)
  
  # Enforce nonnegativity of the infection-age coefficients.
  lower <- rep(0, J)
  upper <- rep(Inf, J)
  
  # Optimize beta with eta fixed using the analytic penalized gradient.
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
  
  # Optionally treat optimizer nonconvergence as an error.
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
  
  # Optionally evaluate the unpenalized NLL for diagnostics.
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
  
  list(
    beta_hat    = fit$par,
    nll_pen     = fit$value,
    nll         = value_nll,
    convergence = fit$convergence,
    message     = fit$message,
    counts      = fit$counts
  )
}


# Alternating estimation --------------------------------------------------


#' Estimate infection-age and healthcare-visit curves
#'
#' Estimates the infection-age and healthcare-visit spline coefficients by
#' alternating between two constrained blockwise optimization steps. Each
#' update minimizes the joint penalized objective with the other coefficient
#' vector held fixed.
#'
#' @param beta_init Numeric vector of length `J` containing the initial
#'   infection-age spline coefficients. Values should be nonnegative.
#' @param eta_init Numeric vector of length `K` containing the initial
#'   healthcare-visit spline coefficients. Values should be nonnegative.
#' @param alpha Nonnegative numeric scalar controlling the common smoothing
#'   penalty applied to \eqn{\beta} and \eqn{\eta}.
#' @param data An `rsv_data` object containing the subject records.
#' @param model An `rsv_model` object containing the daily spline bases,
#'   seasonal RSV circulation curve, and penalty matrices `M_beta` and
#'   `M_eta`.
#' @param rsv_control An `rsv_control` object controlling numerical checks,
#'   clipping, and stabilization.
#' @param subj_pre_list List of subject-level precomputations, with one element
#'   per subject.
#' @param max_outer Positive integer scalar giving the maximum number of
#'   alternating optimization iterations.
#' @param convergence Character scalar specifying the outer convergence
#'   criterion. Supported values are `"objective"`, `"l2"`, and `"both"`.
#' @param tol_obj Nonnegative numeric scalar giving the relative tolerance for
#'   objective-based convergence.
#' @param tol_l2 Nonnegative numeric scalar giving the relative tolerance for
#'   coefficient-based convergence.
#' @param verbose Logical scalar indicating whether to print iteration-level
#'   convergence diagnostics.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{beta_hat}}{Numeric vector of length `J` containing the
#'     estimated infection-age spline coefficients.}
#'   \item{\code{eta_hat}}{Numeric vector of length `K` containing the
#'     estimated healthcare-visit spline coefficients.}
#'   \item{\code{converged}}{Logical scalar indicating whether the selected
#'     outer convergence criterion was satisfied.}
#'   \item{\code{n_outer}}{Integer scalar giving the number of outer
#'     iterations performed.}
#'   \item{\code{inner_conv_beta}}{Logical vector of length `n_outer`
#'     indicating whether each beta-update optimization converged.}
#'   \item{\code{inner_conv_eta}}{Logical vector of length `n_outer`
#'     indicating whether each eta-update optimization converged.}
#'   \item{\code{inner_failed}}{Logical scalar indicating whether any
#'     blockwise optimization reported nonconvergence.}
#'   \item{\code{rel_change_beta}}{Numeric scalar giving the final relative
#'     L2 change in the infection-age coefficients.}
#'   \item{\code{rel_change_eta}}{Numeric scalar giving the final relative
#'     L2 change in the healthcare-visit coefficients.}
#'   \item{\code{rel_change_obj}}{Numeric scalar giving the final relative
#'     change in the joint penalized objective.}
#'   \item{\code{final_objective}}{Numeric scalar giving the joint penalized
#'     objective at the final retained coefficient estimates.}
#' }
#'
#' @details
#' Each outer iteration first updates \eqn{\beta} with \eqn{\eta} fixed, then
#' updates \eqn{\eta} using the new value of \eqn{\beta}. Both blockwise
#' optimizations enforce nonnegative spline coefficients.
#'
#' The joint objective is
#' \deqn{
#'   F(\beta,\eta)
#'   =
#'   \mathrm{NLL}(\beta,\eta)
#'   +
#'   \frac{\alpha}{2}
#'   \left(
#'     \beta^\top M_\beta \beta
#'     +
#'     \eta^\top M_\eta \eta
#'   \right).
#' }
#'
#' With `convergence = "objective"`, the routine stops when the relative
#' objective change is below `tol_obj`. With `"l2"`, both coefficient vectors
#' must have relative L2 changes below `tol_l2`. With `"both"`, both conditions
#' must be satisfied.
#'
#' Inner optimizer convergence is recorded separately from the outer stopping
#' criterion. If a non-finite objective is encountered, alternating estimation
#' stops early and returns the current retained estimates with
#' `converged = FALSE`.
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
  
  # Protect relative-change calculations when a baseline is near zero.
  eps <- 1e-12
  
  beta_curr <- beta_init
  eta_curr  <- eta_init
  
  # Evaluate the joint objective before alternating updates.
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
    
    # Update beta with eta fixed.
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
    
    # Update eta using the new beta estimate (fixed).
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
    
    # Evaluate the joint objective after both block updates.
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
    
    # Compute relative L2 changes using protected baseline norms.
    denom_beta <- sqrt(sum(beta_curr^2))
    if (denom_beta < eps) denom_beta <- 1
    
    denom_eta <- sqrt(sum(eta_curr^2))
    if (denom_eta < eps) denom_eta <- 1
    
    rel_change_beta <- sqrt(sum((beta_new - beta_curr)^2)) / denom_beta
    rel_change_eta  <- sqrt(sum((eta_new  - eta_curr)^2))  / denom_eta
    
    # Compute the relative objective change using a protected baseline.
    denom_obj <- abs(obj_old)
    if (denom_obj < eps) denom_obj <- 1
    
    rel_change_obj <- abs(obj_new - obj_old) / denom_obj
    
    if (verbose) {
      message(sprintf(
        "  rel_obj = %.3e | rel_beta = %.3e | rel_eta = %.3e",
        rel_change_obj, rel_change_beta, rel_change_eta
      ))
    }
    
    beta_curr <- beta_new
    eta_curr  <- eta_new
    obj_old   <- obj_new
    
    # Apply the selected outer convergence criterion.
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
  
  # Retain diagnostics only for completed outer iterations.
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


# Initialization -----------------------------------------------------------


#' Create perturbed initial coefficient values
#'
#' Creates starting values for alternating estimation by applying independent
#' multiplicative perturbations to the infection-age and healthcare-visit
#' spline coefficients. Each perturbed vector is rescaled to preserve the
#' mean coefficient value of the corresponding input vector. When a seed is
#' supplied, the perturbations are reproducible.
#'
#' @param beta Numeric vector of length `J` containing the nonnegative
#'   infection-age spline coefficients to perturb.
#' @param eta Numeric vector of length `K` containing the nonnegative
#'   healthcare-visit spline coefficients to perturb.
#' @param noise_sd Nonnegative numeric scalar giving the standard deviation of
#'   the normal perturbations applied on the log scale.
#' @param seed Optional integer scalar used to initialize the random-number
#'   generator. If supplied, the function calls `set.seed()` and does not
#'   restore the previous RNG state.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{beta_init}}{Numeric vector of length `J` containing the
#'     perturbed infection-age spline coefficients.}
#'   \item{\code{eta_init}}{Numeric vector of length `K` containing the
#'     perturbed healthcare-visit spline coefficients.}
#' }
#'
#' @details
#' For each coefficient vector \eqn{\theta}, the initial values are formed as
#' \deqn{
#'   \theta_j^{\mathrm{init}}
#'   =
#'   \theta_j \exp(\epsilon_j),
#' }
#' where the \eqn{\epsilon_j} are independent normal random variables with
#' mean zero and standard deviation `noise_sd`.
#'
#' Perturbed values below `1e-8` are first bounded at `1e-8`. The resulting
#' vector is then rescaled so that its mean matches the mean of the original
#' coefficient vector. This provides nonnegative starting values near the
#' simulation coefficients without initializing the estimator at the exact
#' values used to generate the data.
#'
#' @examples
#' initial_values <- make_perturbed_initial_values(
#'   beta = c(0.1, 0.2, 0.3),
#'   eta = c(0.3, 0.2, 0.1),
#'   noise_sd = 0.05,
#'   seed = 12345
#' )
#'
#' initial_values$beta_init
#' initial_values$eta_init
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
