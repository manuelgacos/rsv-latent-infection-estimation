#' Subject-level precomputations for RSV model
#'
#' For a given subject \code{subject_i}, this function computes the quantities
#' that depend on the subject and the global model, but do NOT depend on the
#' parameter vectors \code{beta} or \code{eta}. These precomputations can be
#' reused across all evaluations of that subject's contribution to the
#' log-likelihood (e.g., for different values of \code{beta}, \code{eta}).
#'
#' Specifically, this function returns:
#' \itemize{
#'   \item \code{lambda_i} – the subject-shifted calendar-time circulation
#'         curve evaluated on the age grid, i.e.
#'         \eqn{\lambda_i[m] = \lambda(m + B_i)} for \eqn{m = 1,\dots,\text{age\_len}},
#'         where \eqn{B_i} is \code{subject_i$birth_index}.
#'
#'   \item \code{V_i} – the J x 366 cumulative kernel matrix for subject i,
#'         with columns
#'         \eqn{V_i[, 1] = v_i(0) = 0} and
#'         \eqn{V_i[, d + 1] = v_i(d)} for \eqn{d = 1,\dots,365}, as computed
#'         by \code{v_i(lambda_i, phi)} using the precomputed \code{phi}
#'         matrix from the model.
#' }
#'
#' This helper does not depend on \code{beta} or \code{eta}; those parameter
#' vectors enter later through the global precomputations (w_vec, c_vec) and
#' the per-subject likelihood branches.
#'
#' @param subject_i A single subject object, typically created by
#'   \code{make_subject()}, containing at least the fields
#'   \code{birth_index}, \code{visit_age}, and \code{weight}.
#'
#' @param model An \code{"rsv_model"} object created by \code{make_rsv_model()},
#'   which must contain numeric components \code{phi} (J x 365) and
#'   \code{lambda_global} (numeric vector of calendar-time circulation values).
#'
#' @param control An \code{"rsv_control"} object created by
#'   \code{make_rsv_control()}. Its fields \code{check_bounds} and
#'   \code{warn_on_clip} are passed to \code{lambda_shift_i()} to control
#'   how out-of-bounds calendar indices are handled.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{lambda_i} – numeric vector of length equal to the age grid
#'           (typically 365), containing the subject-shifted circulation curve.
#'     \item \code{V_i}      – numeric J x 366 matrix of cumulative kernels for
#'           this subject, as returned by \code{v_i(lambda_i, model$phi)}.
#'   }
#'
#' @export
rsv_precompute_subject <- function(subject_i,
                                   model,
                                   control = make_rsv_control()) {
  # ---- Basic checks on inputs ----------------------------------------------
  if (!is.list(subject_i)) {
    stop("`subject_i` must be a list (typically from make_subject()).")
  }
  if (!inherits(model, "rsv_model")) {
    stop("`model` must be an object of class 'rsv_model' (from make_rsv_model()).")
  }
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("`control` must be an 'rsv_control' object (from make_rsv_control()).")
  }
  
  # Ensure required model components exist
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
  
  # ---- Subject-specific shifted lambda -------------------------------------
  # Uses existing helper: lambda_shift_i(subject_i, model, control)
  lambda_i <- lambda_shift_i(subject_i, model, control)
  
  # Check length consistency with phi's age grid
  if (!is.numeric(lambda_i)) {
    stop("lambda_shift_i() must return a numeric vector.")
  }
  if (length(lambda_i) != ncol(phi)) {
    stop(sprintf(
      "Length of lambda_i (%d) must match ncol(model$phi) (%d).",
      length(lambda_i), ncol(phi)
    ))
  }
  
  # ---- Subject-specific cumulative kernel V_i ------------------------------
  # Uses existing helper: v_i(lambda_i, phi)
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
  
  # ---- Return subject-level precomputations --------------------------------
  list(
    lambda_i = lambda_i,
    V_i      = V_i
  )
}


#' Subject-specific calendar-time shift of the circulation curve
#'
#' @description
#' Returns the subject-specific circulation vector \eqn{\lambda_i(a)} on the
#' *age* grid \eqn{a = 1,\dots,A}, where \eqn{A = ncol(model$B_day)} is the
#' modeled age-window length (typically 365 days).
#'
#' For a subject with birth index \eqn{B_i}, the calendar day corresponding to
#' age \eqn{a} is \eqn{t = B_i + a}. This function therefore maps the age grid
#' to calendar indices \code{cal_idx = B_i + (1:A)} and returns
#' \code{model$lambda_global[cal_idx]}.
#'
#' @details
#' **Indexing convention (important):**
#' \itemize{
#'   \item \code{subject_i$birth_index = B_i} is a 1-based calendar index.
#'   \item \code{age_grid = 1:A} represents ages (in days) from day 1 through day A.
#'   \item The corresponding calendar indices are \code{B_i + 1, ..., B_i + A}.
#' }
#'
#' **Bounds behavior:**
#' \itemize{
#'   \item If \code{control$check_bounds = TRUE}, the function errors if any
#'   computed calendar index is outside \code{1:length(lambda_global)}.
#'   \item Otherwise, the function clips indices to the valid range, optionally
#'   warning if \code{control$warn_on_clip = TRUE}.
#' }
#'
#' @param subject_i A single subject object (from \code{make_subject()}) with
#'   at least \code{birth_index} (integer, 1-based calendar index).
#'
#' @param model An \code{"rsv_model"} object containing:
#'   \itemize{
#'     \item \code{lambda_global}: numeric vector giving \eqn{\lambda(t)} on the
#'           calendar grid, and
#'     \item \code{B_day}: a matrix whose number of columns defines the age
#'           window length \eqn{A} used here.
#'   }
#'
#' @param control An \code{"rsv_control"} object (or compatible list) containing
#'   at least \code{check_bounds} and \code{warn_on_clip}.
#'
#' @return A numeric vector of length \eqn{A = ncol(model$B_day)} giving the
#'   subject-shifted circulation values \eqn{\lambda_i(1),\dots,\lambda_i(A)}.
#'
#' @examples
#' \dontrun{
#' # If birth_index = 1 and A = 365, calendar indices are 2:366
#' lambda_i <- lambda_shift_i(subject_i, model, control)
#' length(lambda_i)  # 365
#' }
#'
#' @export
lambda_shift_i <- function(subject_i, model, control) {
  
  # Extract birth index and global lambda
  B_i <- subject_i$birth_index
  lambda_global <- model$lambda_global
  
  # Age grid length: inferred from B_day columns (assumed already validated)
  age_len  <- ncol(model$B_day)
  age_grid <- seq_len(age_len)
  
  # Calendar indices corresponding to ages 1:age_len
  cal_idx <- B_i + age_grid
  
  # Bounds check only if requested
  if (control$check_bounds) {
    L <- length(lambda_global)
    if (any(cal_idx < 1L | cal_idx > L)) {
      stop(sprintf(
        "lambda_shift_i: calendar indices out of bounds for subject with birth_index = %d.",
        B_i
      ))
    }
  } else {
    # Lenient mode: clip with optional warning
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


#' Subject cumulative accumulator v_i(d)
#'
#' Computes the J x 366 cumulative matrix \eqn{V} for a single subject i, where
#' \eqn{V[, 0] = 0} and \eqn{V[, d] = \sum_{m = 1}^{d} \lambda_{\text{shift}, i}[m] \cdot \phi[, m]}
#' for days \eqn{d = 1, \ldots, 365}. This matches Eqs. (vi_def, vi_element).
#'
#' @param lambda_shift_i Numeric length-365 vector with subject-shifted day
#'   intensities \eqn{\lambda(m + B_i)} for \eqn{m = 1, \ldots, 365}.
#' @param phi Numeric J x 365 matrix whose m-th column is \eqn{\phi(m)} (already
#'   precomputed from your stepwise B-spline basis with eval_rule = "integral").
#'
#' @return Numeric J x 366 matrix \code{V}, where column 1 is \eqn{d = 0} (all zeros),
#'   and column \code{d+1} equals \eqn{v_i(d)} for day \eqn{d = 1, \ldots, 365}.
#'
#' @details
#' Fast path:
#' \itemize{
#'   \item Form \eqn{A_i = \phi \odot \lambda_{\text{shift}, i}} by column scaling:
#'     \code{A <- sweep(phi, 2L, lambda_shift_i, "*")}.
#'   \item Cumulative sum along columns to get \eqn{V[, d]} in a single pass.
#' }
#'
#' Validation identity (useful in tests):
#' \deqn{V[, d] - V[, d-1] = \lambda_{\text{shift}, i}[d] \cdot \phi[, d].}
#'
#' @export
v_i <- function(lambda_shift_i, phi) {
  # --- Shape checks ---
  if (!is.numeric(lambda_shift_i) || length(lambda_shift_i) != 365L) {
    stop("lambda_shift_i must be a numeric vector of length 365.")
  }
  if (!is.matrix(phi) || ncol(phi) != 365L) {
    stop("phi must be a numeric J x 365 matrix (columns are phi[, m]).")
  }
  J <- nrow(phi)
  
  # --- Column scaling: A_i(:, m) = lambda_shift_i[m] * phi(:, m) ---
  A <- sweep(phi, 2L, lambda_shift_i, "*")  # J x 365
  
  # --- Cumulative columns: V[, 0] = 0, V[, d] = V[, d-1] + A[, d] ---
  V <- matrix(0.0, nrow = J, ncol = 366L)
  for (m in 1:365) {
    V[, m + 1L] <- V[, m] + A[, m]
  }
  
  # Optional: dimnames for clarity
  rownames(V) <- rownames(phi)
  colnames(V) <- c("d000", sprintf("d%03d", 1:365))
  V
}


#' Global RSV precomputations for a given (beta, eta)
#'
#' Computes quantities that depend on the parameter vectors (beta, eta)
#' and the global model, but do NOT depend on individual subjects.
#' These can be reused across all subjects in the likelihood.
#'
#' Specifically, this function computes the daily weights
#'   w(m; beta) = beta^T b(m)
#'   c(m; eta)  = eta^T s(m)
#' for m = 1, ..., 365, where the basis matrices B_day and S_day are
#' provided in the rsv_model object. The resulting vectors are then
#' enforced to be nonnegative via the internal helper .enforce_nonneg(),
#' and stabilized for use inside log-likelihood terms via
#' .stabilize_for_log(), which applies a positive floor eps_log from
#' the control object.
#'
#' @param beta Numeric vector of length J (infection-age coefficients).
#'             Must satisfy length(beta) == nrow(model$B_day).
#' @param eta  Numeric vector of length K (detection coefficients).
#'             Must satisfy length(eta) == nrow(model$S_day).
#' @param model An "rsv_model" object created by make_rsv_model(), which
#'   must contain numeric components B_day (J x 365) and S_day (K x 365).
#' @param control An "rsv_control" object created by make_rsv_control().
#'   Its fields tol_clip, check_bounds, warn_on_clip, and eps_log govern how
#'   nonnegativity and log-safety are enforced on w_vec and c_vec.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{w_vec} – numeric length-365 vector w(m; beta), m=1..365,
#'           enforced to be nonnegative and bounded below by eps_log.
#'     \item \code{c_vec} – numeric length-365 vector c(m; eta), m=1..365,
#'           enforced to be nonnegative and bounded below by eps_log.
#'   }
#'
#' @export
rsv_precompute_global <- function(beta,
                                  eta,
                                  model,
                                  control = make_rsv_control()) {
  # ---- Basic checks on model and control ----
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
  
  # Optional: enforce 365-day grid consistency
  if (ncol(B_day) != 365L) {
    stop(sprintf("model$B_day must have 365 columns (found %d).", ncol(B_day)))
  }
  if (ncol(S_day) != 365L) {
    stop(sprintf("model$S_day must have 365 columns (found %d).", ncol(S_day)))
  }
  
  # ---- Dimension consistency with parameters ----
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
  
  # ---- Global per-parameter quantities ----
  w_vec_raw <- w_day(beta, B_day)  # length 365
  c_vec_raw <- c_day(eta,  S_day)  # length 365
  
  # ---- Enforce nonnegativity using unified helper ----
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
  
  # ---- Stabilize for log-safety (apply positive floor eps_log) ----
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
  
  # ---- Return precomputed global terms ----
  list(
    w_vec = w_vec,
    c_vec = c_vec
  )
}


#' Daily hazard weight w(m; beta)
#'
#' Computes the daily hazard weights \eqn{w(m; \beta) = \beta^\top b(m)} 
#' for m = 1, ..., 365. This is the first beta–basis mixer in the RSV model.
#'
#' @param beta Numeric vector of length J. Coefficients for the day basis.
#' @param B_day Numeric matrix of dimension J x 365, whose columns are the
#'   day-basis vectors \eqn{b(m)}. Must match the length of `beta`.
#'
#' @return Numeric vector of length 365 with elements \eqn{w[m] = \beta^\top b(m)}.
#' @examples
#' J <- 4
#' B_day <- matrix(runif(J * 365), nrow = J)
#' beta  <- runif(J)
#' w <- w_day(beta, B_day)
#' stopifnot(length(w) == 365)
#' # Sanity check: manual column sums
#' stopifnot(all.equal(w, colSums(B_day * beta)))
#' @export
w_day <- function(beta, B_day) {
  if (!is.numeric(beta) || !is.numeric(B_day)) {
    stop("beta and B_day must be numeric.")
  }
  if (length(beta) != nrow(B_day)) {
    stop("Length of beta must match number of rows in B_day.")
  }
  # Compute w = B_day^T * beta efficiently (BLAS)
  as.numeric(crossprod(B_day, beta))
}


#' Daily correction weight c(m; eta)
#'
#' Computes the daily weights \eqn{c(m; \eta) = \eta^\top s(m)}
#' for m = 1, ..., 365. This is the second eta–basis mixer in the RSV model.
#'
#' @param eta Numeric vector of length K. Coefficients for the day basis s(m).
#' @param S_day Numeric matrix of dimension K x 365, whose columns are the
#'   day-basis vectors \eqn{s(m)}. Must match the length of `eta`.
#'
#' @return Numeric vector of length 365 with elements \eqn{c[m] = \eta^\top s(m)}.
#' @examples
#' K <- 4
#' S_day <- matrix(runif(K * 365), nrow = K)
#' eta   <- runif(K)
#' c_vec <- c_day(eta, S_day)
#' stopifnot(length(c_vec) == 365)
#' # Sanity check: manual column sums
#' stopifnot(all.equal(c_vec, colSums(S_day * eta)))
#' @export
c_day <- function(eta, S_day) {
  if (!is.numeric(eta) || !is.numeric(S_day)) {
    stop("eta and S_day must be numeric.")
  }
  if (length(eta) != nrow(S_day)) {
    stop("Length of eta must match number of rows in S_day.")
  }
  # Compute c = S_day^T * eta efficiently (BLAS)
  as.numeric(crossprod(S_day, eta))
}


#' Internal helper: enforce nonnegativity on weight-like vectors
#'
#' Ensures that a numeric vector \code{x} satisfies x >= 0 up to a tolerance.
#' Values in [-tol_clip, 0) are treated as numerical noise and clipped to 0.
#' Values < -tol_clip are treated as serious violations:
#'   * if check_bounds = TRUE: an error is thrown;
#'   * if check_bounds = FALSE: they are clipped to 0 with a warning.
#'
#' This is intended for weight-like quantities such as w_vec and c_vec that
#' are constrained by the model to be nonnegative, but are not probabilities
#' (so we do not impose an upper bound here).
#'
#' @param x Numeric vector to be checked and corrected.
#' @param name Character string used in warnings/errors (e.g. "w_vec", "c_vec").
#' @param control An rsv_control object from make_rsv_control(), providing
#'   tol_clip, check_bounds, and warn_on_clip.
#'
#' @return A numeric vector with the same length as x, with all negative
#'   entries replaced by 0 (subject to error behavior when check_bounds = TRUE
#'   and serious violations are present).
#'
#' @keywords internal
.enforce_nonneg <- function(x,
                            name    = "value",
                            control = make_rsv_control()) {
  # ---- Basic checks ----
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric in .enforce_nonneg().", name))
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` contains non-finite values (NA, NaN, or Inf).", name))
  }
  
  # Extract control fields with simple fallbacks
  tol_clip     <- if (!is.null(control$tol_clip))     control$tol_clip     else 1e-12
  check_bounds <- if (!is.null(control$check_bounds)) control$check_bounds else FALSE
  warn_on_clip <- if (!is.null(control$warn_on_clip)) control$warn_on_clip else TRUE
  
  if (!is.numeric(tol_clip) || length(tol_clip) != 1L || tol_clip < 0 || !is.finite(tol_clip)) {
    stop("control$tol_clip must be a single nonnegative finite numeric value.")
  }
  
  # ---- Identify negatives ----
  neg_idx <- which(x < 0)
  if (length(neg_idx) == 0L) {
    # Fast path: nothing to do
    return(x)
  }
  
  x_neg <- x[neg_idx]
  
  # Tiny negatives: within [-tol_clip, 0)
  tiny_idx    <- neg_idx[x_neg >= -tol_clip]
  # Serious negatives: < -tol_clip
  serious_idx <- neg_idx[x_neg < -tol_clip]
  
  # ---- Serious violations: possibly error ----
  if (length(serious_idx) > 0L && isTRUE(check_bounds)) {
    min_val <- min(x[serious_idx])
    stop(sprintf(
      "%s has %d values < -tol_clip (min = %g) violating nonnegativity.",
      name, length(serious_idx), min_val
    ))
  }
  
  # ---- Perform clipping ----
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
    # For serious violations, always warn in lenient mode
    min_val <- min(x_neg[x_neg < -tol_clip])
    warning(sprintf(
      "%s had %d values < -tol_clip (min = %g); clipped to 0. This may indicate an issue with the model or parameters.",
      name, length(serious_idx), min_val
    ))
  }
  
  x
}


#' Internal helper: stabilize nonnegative values for use inside log()
#'
#' Given a numeric vector \code{x} that is expected to be nonnegative, this
#' helper enforces a strictly positive lower bound \code{eps_log} so that
#' \code{log(x)} is numerically well-defined and finite. Values in
#' \code{[0, eps_log)} are "floored" to \code{eps_log}. Negative values are
#' treated as violations when \code{check_bounds = TRUE}.
#'
#' This helper is intended for quantities like the day-wise weights
#' \code{w(m; beta)} or \code{c(m; eta)} that are nonnegative by construction
#' but may become very small or slightly negative due to numerical error.
#'
#' @param x Numeric vector, expected to be nonnegative.
#' @param name Character string used in warning/error messages (e.g.,
#'   \code{"w_vec"}, \code{"c_vec"}).
#' @param eps_log Positive numeric scalar. Strictly positive floor used to
#'   stabilize \code{x} for use inside \code{log()}. Typically taken from
#'   \code{control$eps_log}.
#' @param check_bounds Logical; if \code{TRUE}, any negative values in \code{x}
#'   are treated as errors. If \code{FALSE}, negative values are still
#'   stabilized but no error is thrown.
#' @param warn_on_clip Logical; if \code{TRUE}, emit a warning whenever one or
#'   more entries of \code{x} are raised to \code{eps_log}.
#'
#' @return A numeric vector of the same length as \code{x}, with all entries
#'   satisfying \code{x >= eps_log}.
#'
#' @keywords internal
.stabilize_for_log <- function(x,
                               name,
                               eps_log,
                               check_bounds = FALSE,
                               warn_on_clip = TRUE) {
  
  # Basic checks
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric in .stabilize_for_log().", name))
  }
  if (!is.numeric(eps_log) || length(eps_log) != 1L ||
      !is.finite(eps_log) || eps_log <= 0) {
    stop("`eps_log` must be a single positive finite numeric value in .stabilize_for_log().")
  }
  
  # Detect negative values
  neg_idx <- which(x < 0)
  if (length(neg_idx) > 0L) {
    if (check_bounds) {
      stop(sprintf(
        "%s has %d negative values in .stabilize_for_log(); input is expected to be nonnegative.",
        name, length(neg_idx)
      ))
    } else {
      # Lenient mode: raise negatives to eps_log
      if (warn_on_clip) {
        warning(sprintf(
          "%s has %d negative values; stabilizing them to eps_log = %g in .stabilize_for_log().",
          name, length(neg_idx), eps_log
        ))
      }
      x[neg_idx] <- eps_log
    }
  }
  
  # Floor small nonnegative values to eps_log
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


#' Subject cumulative hazard H_i(m; beta)
#'
#' Computes the subject- and day-specific cumulative hazard
#' \deqn{H_i(m; \beta) = \beta^\top v_i(m)} for m = 1, ..., 365,
#' where columns of V_i are the cumulative kernels v_i(m) (Card #2).
#'
#' @param beta Numeric vector of length J. Coefficients for the day basis.
#' @param V_i Numeric matrix J x 366. Column d is v_i(d), with column 1
#'   representing v_i(0) = 0 (by convention). If your V_i has 366 columns,
#'   this function uses columns 2:366 to return H_i(1:365).
#' @param include_day0 Logical, default FALSE. If TRUE, prepend H_i(0)=0.
#'
#' @return Numeric vector of length 365 (or 366 if include_day0=TRUE).
#' @examples
#' J <- 4
#' # Build a fake cumulative V_i: v_i(0)=0, then cumulative sums of random Jx365
#' A  <- matrix(rexp(J * 365, rate = 1), nrow = J)
#' V  <- cbind(0, t(apply(A, 1, cumsum)))  # WRONG shape; fix to J x 366
#' V_i <- matrix(0, nrow = J, ncol = 366)
#' V_i[, 1] <- 0
#' V_i[, 2:366] <- apply(A, 1, cumsum)
#' beta <- runif(J)
#' H <- H_i(beta, V_i)
#' stopifnot(length(H) == 365)
#' # Manual check for a single day m:
#' m <- 123
#' stopifnot(all.equal(H[m], sum(beta * V_i[, m + 1])))
#' @export
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
  # Use columns 2:366 → v_i(1:365)
  H <- as.numeric(crossprod(V_i[, 2:366, drop = FALSE], beta))
  if (isTRUE(include_day0)) {
    H <- c(0, H)
  }
  H
}


#' Subject survival curve \bar F_i(m; beta)
#'
#' Computes the subject-specific survival values:
#'   Fbar_i(m) = exp( - H_i(m) )
#' for m = 1,...,365 (or including m = 0 if include_day0 = TRUE).
#'
#' Numerical clipping is applied to ensure that survival values remain within
#' (0, 1] up to floating-point tolerance. This protects against underflow
#' (exp(-H) = 0) and slight numerical drift above 1.
#'
#' @param H Numeric vector of length 365 (if include_day0 = FALSE)
#'   or length 366 (if include_day0 = TRUE). Must contain the cumulative
#'   hazard values H_i(m).
#' @param include_day0 Logical, default FALSE. If TRUE, expects H[1] = H_i(0)=0
#'   and returns the corresponding survival value Fbar_i(0)=1.
#' @param check_bounds Logical; if TRUE, treat severe violations of the
#'   theoretical bounds (values far outside (0,1]) as errors. Intended for
#'   debugging and model validation.
#' @param warn_on_clip Logical; if TRUE, emit a warning whenever clipping is
#'   applied (useful for diagnosing small numerical drift).
#'
#' @return Numeric vector of survival values with same length as H.
#'
#' @export
Fbar_from_H <- function(H,
                        include_day0 = FALSE,
                        check_bounds = FALSE,
                        warn_on_clip = TRUE) {
  
  if (!is.numeric(H))
    stop("H must be numeric.")
  
  len_H <- length(H)
  
  if (include_day0) {
    # Expect H[1] = H_i(0)
    if (len_H != 366)
      stop("If include_day0=TRUE, H must have length 366 (H_i(0)...H_i(365)).")
    
    Fbar <- exp(-H)
    
  } else {
    # Expect length 365
    if (len_H != 365)
      stop("If include_day0=FALSE, H must have length 365 (H_i(1)...H_i(365)).")
    
    Fbar <- exp(-H)
  }
  
  # --- Numerical clipping for stability ---
  # Survival values must lie in (0,1].
  # Use a small positive lower bound to avoid underflow to 0.
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


#' Convenience wrapper for \bar F_i(m; beta)
#'
#' Computes \eqn{\bar F_i(m; \beta)} by first computing H_i(m; beta) from
#' `beta` and `V_i`, then applying exp(-H).
#'
#' @param beta Numeric vector length J.
#' @param V_i  Numeric matrix J x 366, with column 1 being v_i(0)=0.
#' @param include_day0 Logical, default FALSE. If TRUE, prepend day 0.
#' @param check_bounds Logical; if TRUE, propagate strict bound checking to
#'   the internal survival computation.
#' @param warn_on_clip Logical; if TRUE, propagate clipping warnings from the
#'   internal survival computation.
#'
#' @return Numeric vector of length 365 (or 366 if include_day0=TRUE).
#' @examples
#' # Build a small synthetic example
#' set.seed(1)
#' J <- 3
#' A  <- matrix(abs(rnorm(J * 365)), nrow = J)
#' V_i <- matrix(0, nrow = J, ncol = 366)
#' V_i[, 2:366] <- apply(A, 1, cumsum)
#' beta <- runif(J)
#' Fbar <- Fbar_i(beta, V_i)
#' stopifnot(length(Fbar) == 365, all(Fbar >= 0), all(Fbar <= 1))
#' @export
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


#' Event probability per day: pi_i(m; beta)
#'
#' Computes the subject- and day-specific event probability
#' \deqn{\pi_i(m; \beta) = 1 - \exp\{-\lambda_i[m] \cdot w[m]\}, \quad m=1,\dots,365,}
#' where \eqn{w[m] = \beta^\top b(m)} and \eqn{\lambda_i[m] = \lambda(m+B_i)}.
#'
#' Prefer calling `pi_from_w(lambda_shift_i, w)` if `w` is already computed
#' (e.g., once per beta using `w_day(beta, B_day)`). The wrapper `pi_i()`
#' will compute `w` internally from `beta` and `B_day` for convenience.
#'
#' Numerical clipping is applied to ensure that the resulting probabilities
#' lie in [0, 1] up to floating-point tolerance. This protects against small
#' numerical drift outside [0, 1].
#'
#' @param lambda_shift_i Numeric vector length 365 with \eqn{\lambda(m+B_i)} for subject i.
#' @param w Numeric vector length 365 with \eqn{w[m] = \beta^\top b(m)}.
#' @param check_bounds Logical; if TRUE, severe violations of the [0,1] bounds
#'   (values far outside [0,1] beyond a tolerance) trigger an error. Intended
#'   for debugging and model validation.
#' @param warn_on_clip Logical; if TRUE, emit a warning whenever clipping is
#'   applied to the probabilities.
#'
#' @return Numeric vector length 365 with \eqn{\pi_i[m]} in [0,1] up to numerical
#'   tolerance.
#' @examples
#' # Core usage with precomputed w:
#' set.seed(1)
#' lambda_i <- rexp(365, 0.2)   # shifted lambda for subject i
#' w <- runif(365)              # w(m; beta)
#' pi <- pi_from_w(lambda_i, w)
#' stopifnot(length(pi) == 365, all(pi >= 0), all(pi <= 1))
#' @export
pi_from_w <- function(lambda_shift_i,
                      w,
                      check_bounds = FALSE,
                      warn_on_clip = TRUE) {
  if (!is.numeric(lambda_shift_i) || length(lambda_shift_i) != 365L)
    stop("lambda_shift_i must be numeric length 365.")
  if (!is.numeric(w) || length(w) != 365L)
    stop("w must be numeric length 365.")
  
  # Stable form: 1 - exp(-x) == -expm1(-x)
  x  <- lambda_shift_i * w
  pi <- -expm1(-x)
  
  # Numerical guard to keep probabilities within [0,1] up to tolerance
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


#' Q_i(m; beta) = \bar F_i(m-1) * pi_i(m)
#'
#' Vectorized computation of the visit-day kernel Q_i over m = 1..365, given
#' survival tail \bar F_i(d; beta) for d = 0..365 and per-day visit prob
#' pi_i(m; beta) for m = 1..365.
#'
#' @param Fbar_i Numeric length-366 vector: \bar F_i(d; beta) for d = 0..365.
#'               The first entry corresponds to d = 0.
#' @param pi_i   Numeric length-365 vector: pi_i(m; beta) for m = 1..365.
#' @param check_bounds Logical; if TRUE, propagate strict bound checking to
#'   the internal probability computations.
#' @param warn_on_clip Logical; if TRUE, emit warnings whenever clipping is
#'   applied to Fbar_i or pi_i.
#'
#' @return Numeric length-365 vector Q_i with Q_i[m] = Fbar_i[m] * pi_i[m],
#'         i.e., \bar F_i(m-1) * pi_i(m).
#' @export
Q_i <- function(Fbar_i,
                pi_i,
                check_bounds = FALSE,
                warn_on_clip = TRUE) {
  # --- shape checks ---
  if (!is.numeric(Fbar_i) || length(Fbar_i) != 366L) {
    stop("Fbar_i must be a numeric vector of length 366 (d = 0..365).")
  }
  if (!is.numeric(pi_i) || length(pi_i) != 365L) {
    stop("pi_i must be a numeric vector of length 365 (m = 1..365).")
  }
  
  # --- stability guards using unified clipping ---
  # Fbar must lie in (0, 1], pi must lie in [0, 1].
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
  
  # --- vectorized product with the (m-1) shift for Fbar ---
  # For m=1..365, use Fbar_i[m] (which is \bar F at d=m-1)
  # i.e., Fbar_i[m] corresponds to \bar F_i(m-1).
  Q <- Fbar_i[1:365] * pi_i
  
  # Final clipping to [0,1] to counter tiny numerical drift
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


#' U_i(beta) = sum_m s(m) * Q_i(m; beta)
#'
#' Computes the K-vector U_i by multiplying the precomputed day-basis S_day
#' (K x 365; columns are s(m)) by the Q_i vector (length 365).
#'
#' @param Q_i   Numeric length-365 vector: Q_i(m; beta) for m = 1..365.
#' @param S_day Numeric K x 365 matrix: columns are s(m).
#' @param check_bounds Logical; if TRUE, apply strict bound checking to Q_i
#'   via the internal clipping helper.
#' @param warn_on_clip Logical; if TRUE, emit a warning whenever clipping is
#'   applied to Q_i.
#'
#' @return Numeric length-K vector U_i = S_day %*% Q_i.
#' @export
U_i <- function(Q_i,
                S_day,
                check_bounds = FALSE,
                warn_on_clip = TRUE) {
  ## --- shape checks ---
  if (!is.numeric(Q_i) || length(Q_i) != 365L) {
    stop("Q_i must be a numeric vector of length 365.")
  }
  
  # Allow base matrix or Matrix package classes
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
  
  ## --- stability guard for Q_i using unified clipping ---
  # Q_i should be in [0, 1] since it is built from probabilities.
  Q_i <- .clip_to_range(
    Q_i,
    lower        = 0.0,
    upper        = 1.0,
    name         = "Q_i",
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  ## --- compute U_i = S_day %*% Q_i ---
  U <- as.vector(S_day %*% Q_i)
  
  ## --- clean tiny negatives from numerical noise ---
  U[U < 0 & U > -1e-12] <- 0
  
  U
}


#' tau_i(beta, eta) = 1 - t(eta) %*% U_i(beta)
#'
#' Computes the scalar \eqn{\tau_i(\beta, \eta) = 1 - \eta^\top U_i(\beta)}.
#' Numerically guards against values at or below a small threshold \code{eps_tau}.
#'
#' @param U_i Numeric vector of length K: U_i(beta).
#' @param eta Numeric vector of length K: visit-parameter vector \eqn{\eta}.
#' @param eps_tau Nonnegative numeric scalar; minimum allowed value for tau.
#'   Default \code{1e-12}. Set to \code{0} to disable clipping.
#' @param check_bounds Logical; if TRUE, throw an error when tau violates
#'   the admissible interval [eps_tau, 1].
#' @param warn_on_clip Logical; if TRUE, emit a warning when tau is clipped
#'   to the admissible interval [eps_tau, 1].
#'
#' @return Numeric scalar tau in \eqn{[eps_\tau, 1]}.
#' @export
tau_i <- function(U_i,
                  eta,
                  eps_tau = 1e-12,
                  check_bounds = FALSE,
                  warn_on_clip = TRUE) {
  # --- shape checks ---
  if (!is.numeric(U_i) || !is.numeric(eta))
    stop("U_i and eta must be numeric.")
  if (length(U_i) != length(eta))
    stop("U_i and eta must have the same length (K).")
  if (!is.finite(eps_tau) || eps_tau < 0)
    stop("eps_tau must be a nonnegative finite number.")
  
  # --- core computation ---
  dot <- sum(eta * U_i)
  if (!is.finite(dot))
    stop("Non-finite dot product: check inputs.")
  
  tau <- 1 - dot
  
  # --- clipping / diagnostics ---
  # Valid range is [eps_tau, 1].
  # tau_i should not use .Machine$double.eps as a lower bound,
  # because extremely small tau leads to numerical instability in
  # normalization and log-likelihood evaluation. A model-appropriate floor
  # (eps_tau, default 1e-12) is used instead.
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


#' Internal helper: clip numeric values into a bounded interval
#'
#' Clips a numeric vector \code{x} into the interval \code{[lower, upper]},
#' with optional strict boundary checking and optional warnings when clipping
#' occurs. This utility is intended for internal use in functions that compute
#' probabilities, survival values, or other quantities that are theoretically
#' bounded within a known interval.
#'
#' The behavior is controlled by two arguments:
#' \itemize{
#'   \item \code{check_bounds}: If \code{TRUE}, values that lie \emph{far}
#'   outside the interval \code{[lower, upper]} (beyond a tolerance \code{tol})
#'   trigger an error. This is useful for debugging or strict validation.
#'
#'   \item \code{warn_on_clip}: If \code{TRUE}, a warning is emitted whenever
#'   clipping occurs, even if values lie only slightly outside the interval
#'   (within \code{tol}). This helps detect small numerical drift while keeping
#'   computation stable.
#' }
#'
#' The final output is always clipped to \code{[lower, upper]} regardless of
#' the settings, unless \code{check_bounds = TRUE} and severe violations are
#' detected, in which case an error is raised.
#'
#' @param x Numeric vector. Values to be clipped into \code{[lower, upper]}.
#' @param lower,upper Numeric scalars defining the clipping interval. Must
#'   satisfy \code{lower <= upper}.
#' @param name Character string used in error and warning messages (e.g.,
#'   \code{"pi_i"}, \code{"Q_i"}, \code{"Fbar_i"}).
#' @param check_bounds Logical; if \code{TRUE}, severe violations (values below
#'   \code{lower - tol} or above \code{upper + tol}) yield an error.
#' @param warn_on_clip Logical; if \code{TRUE}, emit a warning whenever any
#'   value is clipped.
#' @param tol Numeric tolerance used to distinguish slight floating-point drift
#'   from genuine violations. Must be nonnegative.
#'
#' @return A numeric vector of the same length as \code{x}, with all values
#'   clipped into \code{[lower, upper]}.
#'
#' @keywords internal
#' @examples
#' x <- c(-1e-14, 0.5, 1 + 2e-14)
#' .clip_to_range(x, 0, 1, name = "example")
#'
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
  
  # Identify any out-of-bounds values
  below <- x < lower
  above <- x > upper
  any_oob <- any(below | above)
  
  if (!any_oob) {
    # Fast path: no clipping needed
    return(x)
  }
  
  # Severe violations beyond numerical tolerance
  far_below <- x < (lower - tol)
  far_above <- x > (upper + tol)
  n_far <- sum(far_below | far_above)
  
  if (check_bounds && n_far > 0L) {
    stop(sprintf(
      "%s has %d values far outside [%g, %g] (beyond tol = %g).",
      name, n_far, lower, upper, tol
    ))
  }
  
  # Warn when clipping occurs
  if (warn_on_clip) {
    n_oob <- sum(below | above)
    warning(sprintf(
      "%s has %d values slightly outside [%g, %g]; clipping applied.",
      name, n_oob, lower, upper
    ))
  }
  
  # Clip and return
  pmin(pmax(x, lower), upper)
}