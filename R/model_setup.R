#' Construct a control object for RSV likelihood computations
#'
#' This function creates a small list of control parameters that govern
#' numerical checks, clipping behavior, and related diagnostics in the RSV
#' likelihood and its internal building blocks. It does not affect the
#' statistical model itself; instead, it controls how strictly intermediate
#' quantities are validated and how numerical edge cases are handled.
#'
#' Typical usage:
#' \itemize{
#'   \item During development or debugging, use the defaults:
#'     \code{make_rsv_control(check_bounds = TRUE, warn_on_clip = TRUE)}.
#'
#'   \item For large-scale simulation or production runs, you may turn off
#'     most checks and warnings for speed and cleaner output:
#'     \code{make_rsv_control(check_bounds = FALSE, warn_on_clip = FALSE)}.
#' }
#'
#' The fields of the returned list are intended to be passed to internal
#' functions (e.g., \code{loglik_i()}, \code{tau_i()}, \code{pi_from_w()})
#' and then forwarded to utilities such as \code{.clip_to_range()}.
#'
#' @param check_bounds Logical. If \code{TRUE}, internal functions are allowed
#'   (and encouraged) to perform strict bound checks on probabilities,
#'   survival values, and related quantities. Severe violations may trigger
#'   errors instead of being silently corrected.
#'
#' @param warn_on_clip Logical. If \code{TRUE}, internal functions emit
#'   warnings whenever clipping is applied to a quantity that should, in
#'   principle, lie within a fixed interval (e.g., probabilities in [0, 1],
#'   or \eqn{\tau_i} being below a minimum threshold). This is particularly
#'   useful for detecting small numerical drift during development.
#'
#' @param eps_tau Numeric scalar. Minimum allowed value for quantities like
#'   \eqn{\tau_i}, where zero or negative values would cause numerical issues
#'   (for example, inside logarithms). Typical default is \code{1e-12}.
#'
#' @param tol_clip Numeric scalar. Default tolerance passed to clipping
#'   helpers (such as \code{.clip_to_range()}) to distinguish floating-point
#'   noise from genuine bound violations. Must be nonnegative.
#'
#' @param eps_log Numeric scalar. Minimum allowed value for nonnegative
#'   quantities that will be used inside logarithms (for example, the
#'   day-wise weights \code{w(m; beta)} or \code{c(m; eta)}). Values in
#'   \code{[0, eps_log)} are stabilized upward to \code{eps_log} to avoid
#'   \code{log(0)} while keeping the effect localized to very small values.
#'
#' @return A list with class \code{"rsv_control"} containing:
#'   \itemize{
#'     \item \code{check_bounds} – logical flag for strict bound checking.
#'     \item \code{warn_on_clip} – logical flag for warnings on clipping.
#'     \item \code{eps_tau}      – minimum value for \eqn{\tau_i}-like terms.
#'     \item \code{tol_clip}     – default tolerance for clipping utilities.
#'     \item \code{eps_log}      – minimum value for log-stabilized
#'                                nonnegative quantities.
#'   }
#'
#' The returned object is designed to be passed unchanged throughout the
#' likelihood computation and should be treated as read-only.
#'
make_rsv_control <- function(
    check_bounds = TRUE,
    warn_on_clip = TRUE,
    eps_tau      = 1e-12,
    tol_clip     = 1e-12,
    eps_log      = 1e-12
) {
  
  # ---- Basic input validation -----------------------------------------------
  
  if (!is.logical(check_bounds) || length(check_bounds) != 1L || is.na(check_bounds)) {
    stop("`check_bounds` must be a single non-NA logical value.")
  }
  
  if (!is.logical(warn_on_clip) || length(warn_on_clip) != 1L || is.na(warn_on_clip)) {
    stop("`warn_on_clip` must be a single non-NA logical value.")
  }
  
  if (!is.numeric(eps_tau) || length(eps_tau) != 1L || !is.finite(eps_tau) || eps_tau <= 0) {
    stop("`eps_tau` must be a single positive finite numeric value.")
  }
  
  if (!is.numeric(tol_clip) || length(tol_clip) != 1L || !is.finite(tol_clip) || tol_clip < 0) {
    stop("`tol_clip` must be a single nonnegative finite numeric value.")
  }
  
  if (!is.numeric(eps_log) || length(eps_log) != 1L || !is.finite(eps_log) || eps_log <= 0) {
    stop("`eps_log` must be a single positive finite numeric value.")
  }
  
  # ---- Construct control object --------------------------------------------
  
  out <- list(
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip,
    eps_tau      = eps_tau,
    tol_clip     = tol_clip,
    eps_log      = eps_log
  )
  
  class(out) <- "rsv_control"
  
  return(out)
}


# NOTE: Was updated in V6
#' RSV simulation model with user-specified calendar-time lambda
#'
#' Builds an \code{rsv_model} for simulations using a user-supplied calendar-time
#' circulation curve \eqn{\lambda(t)} over a calendar grid of length \code{calendar_len}.
#'
#' @param lambda_global Numeric vector specifying the calendar-time lambda curve.
#'   Must be nonnegative and of length \code{calendar_len}.
#' @param calendar_len Integer. Total number of calendar days. Must equal
#'   \code{length(lambda_global)}.
#' @param J Integer. Number of basis functions for the infection-age curve \eqn{w}.
#' @param K Integer. Number of basis functions for the visit-age curve \eqn{c}.
#' @param degree Integer. Spline degree (typically 3 for cubic).
#' @param days Integer. Number of days in the age grid used to build \code{B_day}, 
#'   \code{S_day}, and \code{phi}. For this simulation setup we enforce \code{days == age_len}.
#' @param age_len Integer. Maximum age allowed in the simulation (e.g., 365 days). Must equal \code{days}.
#' @param birth_max Integer. Maximum birth day index used when sampling birthdays.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{model} – the \code{rsv_model} object (via \code{make_rsv_model()})
#'     \item \code{B_day}, \code{S_day}, \code{phi} – basis matrices
#'     \item \code{lambda_global}, \code{calendar_len}, \code{birth_max}, \code{age_len}
#'     \item \code{lambda_params} – a named list describing the lambda input
#'   }
#'
#' @export
make_sim_model_lambda_custom <- function(lambda_global,
                                         J = 18L,
                                         K = 18L,
                                         degree = 3L,
                                         days = 365L,
                                         age_len = 365L,
                                         birth_max = 365L,
                                         calendar_len = 730L) {
  # ---- Basic validation ----
  stopifnot(is.numeric(lambda_global))
  stopifnot(is.numeric(calendar_len), length(calendar_len) == 1L)
  stopifnot(calendar_len > 0, length(lambda_global) == calendar_len)
  stopifnot(all(lambda_global >= 0))
  
  stopifnot(J > 0L, K > 0L, days > 0L, age_len > 0L, birth_max > 0L)
  
  if (as.integer(days) != as.integer(age_len)) {
    stop(sprintf("make_sim_model_lambda_custom: require days == age_len; got days=%d, age_len=%d.",
                 as.integer(days), as.integer(age_len)))
  }
  
  required_calendar <- as.integer(birth_max) + as.integer(age_len)
  if (calendar_len < required_calendar) {
    stop(sprintf(
      paste0("make_sim_model_lambda_custom: calendar_len too short. ",
             "Need calendar_len >= birth_max + age_len = %d, got %d."),
      required_calendar, calendar_len
    ))
  }
  
  # ---- Basis construction ----
  bases <- make_sim_bases(J = J, K = K, degree = degree, days = as.integer(age_len))
  B_day <- bases$B_day
  S_day <- bases$S_day
  phi   <- bases$phi
  
  # ---- Model object ----
  model <- make_rsv_model(
    B_day         = B_day,
    S_day         = S_day,
    phi           = phi,
    lambda_global = lambda_global
  )
  
  model$D_beta <- make_D2(J)
  model$M_beta <- crossprod(model$D_beta)
  
  model$D_eta  <- make_D2(K)
  model$M_eta  <- crossprod(model$D_eta)
  
  # ---- Return ----
  list(
    model         = model,
    B_day         = B_day,
    S_day         = S_day,
    phi           = phi,
    lambda_global = lambda_global,
    calendar_len  = as.integer(calendar_len),
    birth_max     = as.integer(birth_max),
    age_len       = as.integer(age_len),
    lambda_params = list(
      type   = "custom",
      source = "user_provided"
    )
  )
}


#' Construct day-level spline bases for RSV simulations
#'
#' @description
#' Builds the day-level basis matrices used throughout the RSV simulation
#' pipeline:
#' \itemize{
#'   \item \code{B_day}: basis for the infection-age density \eqn{w(m;\beta)}
#'   \item \code{S_day}: basis for the visit/detection function \eqn{c(m;\eta)}
#'   \item \code{phi}:  day-integrated basis used to construct cumulative kernels
#'         \eqn{v_i(d)} (via \code{v_i(lambda_i, phi)}) in subject precomputations
#' }
#'
#' The basis is constructed on the discrete day grid \eqn{m = 1,\dots,\code{days}}
#' using B-splines of degree \code{degree} with boundary knots
#' \code{c(0, days)}. Interior knots are created via \code{make_interior_knots()},
#' using either quantile-based or equally spaced placement over the day grid.
#'
#' @details
#' The infection-age basis (\code{B_day}) is built using interior knots for a
#' \eqn{J}-dimensional spline basis. The same basis is also used to build
#' \code{phi}, which represents day-integrated basis columns
#' \eqn{\phi(m) = \int_{m-1}^{m} b(x)\,dx}. If \code{K == J}, the visit-age basis
#' \code{S_day} is set equal to \code{B_day}; otherwise, a separate \eqn{K}-dimensional
#' spline basis is built for \code{S_day}.
#'
#' This function is meant for simulation setup (Bucket A) and returns plain
#' matrices plus metadata needed to construct an \code{rsv_model} object.
#'
#' @param J Integer > 0. Number of basis functions for the infection-age density
#'   \eqn{w(m;\beta)}. Determines \code{nrow(B_day)} and \code{nrow(phi)}.
#' @param K Integer > 0. Number of basis functions for the visit/detection
#'   function \eqn{c(m;\eta)}. Determines \code{nrow(S_day)}.
#' @param degree Integer >= 0. Degree of the B-spline basis (e.g., \code{3L} for cubic).
#' @param days Integer > 0. Number of days in the age grid. Determines the number
#'   of columns in \code{B_day}, \code{S_day}, and \code{phi}.
#' @param knot_method Character string. Method for placing interior knots:
#'   \code{"quantile"} (default) or \code{"equispaced"}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{B_day}}{Numeric matrix of dimension \code{J x days}. Column \code{m}
#'     is the day-basis vector \eqn{b(m)} used to form \eqn{w(m;\beta)=\beta^\top b(m)}.}
#'   \item{\code{S_day}}{Numeric matrix of dimension \code{K x days}. Column \code{m}
#'     is the day-basis vector \eqn{s(m)} used to form \eqn{c(m;\eta)=\eta^\top s(m)}.}
#'   \item{\code{phi}}{Numeric matrix of dimension \code{J x days}. Column \code{m}
#'     represents the day-integrated basis \eqn{\phi(m)=\int_{m-1}^{m} b(x)\,dx},
#'     used in cumulative kernel construction.}
#'   \item{\code{Boundary.knots}}{Numeric length-2 vector \code{c(0, days)} giving the
#'     spline domain endpoints.}
#'   \item{\code{degree}}{Integer spline degree used.}
#'   \item{\code{days}}{Integer number of days used.}
#' }
#'
#' @seealso \code{\link{make_interior_knots}}, \code{\link{build_B_day_from_spline}},
#'   \code{\link{make_rsv_model}}
#'
#' @export
make_sim_bases <- function(J = 18L,
                           K = 18L,
                           degree = 3L,
                           days = 365L,
                           knot_method = c("quantile", "equispaced")) {
  knot_method <- match.arg(knot_method)
  stopifnot(J > 0L, K > 0L, days > 0L)
  
  day_grid <- seq_len(days)
  Boundary.knots <- c(0, days)  # domain for the spline basis
  
  # ---- Interior knots for w (infection-age basis) ----
  knots_w <- make_interior_knots(
    J              = J,
    degree         = degree,
    Boundary.knots = Boundary.knots,
    method         = knot_method,
    x_grid         = day_grid
  )
  
  # ---- B_day: stepwise basis for w (used in w_day) ----
  # Columns: b(m), m = 1..365; dimension J x 365
  B_day <- build_B_day_from_spline(
    knots          = knots_w,
    degree         = degree,
    J              = J,
    eval_rule      = "integral",   # evaluate at midpoints m - 0.5
    Boundary.knots = Boundary.knots
  )
  
  # ---- phi: integrated basis for v_i (used in cumulative kernel) ----
  # Columns: phi(m) = ∫_{m-1}^m b(x) dx, m = 1..365; dimension J x 365
  phi <- build_B_day_from_spline(
    knots          = knots_w,
    degree         = degree,
    J              = J,
    eval_rule      = "integral",   # integrated B-splines → day integrals
    Boundary.knots = Boundary.knots
  )
  
  # ---- S_day: basis for c (visit-age basis) ----
  # For now we mirror the w-basis; if K != J we rebuild with K.
  if (K == J) {
    S_day <- B_day
  } else {
    knots_c <- make_interior_knots(
      J              = K,
      degree         = degree,
      Boundary.knots = Boundary.knots,
      method         = knot_method,
      x_grid         = day_grid
    )
    S_day <- build_B_day_from_spline(
      knots          = knots_c,
      degree         = degree,
      J              = K,
      eval_rule      = "integral",
      Boundary.knots = Boundary.knots
    )
  }
  
  # ---- Basic sanity checks ----
  stopifnot(
    nrow(B_day) == J, ncol(B_day) == days,
    nrow(S_day) == K, ncol(S_day) == days,
    nrow(phi)   == J, ncol(phi)   == days
  )
  
  list(
    B_day          = B_day,
    S_day          = S_day,
    phi            = phi,
    Boundary.knots = Boundary.knots,
    degree         = degree,
    days           = days
  )
}


#' Construct interior knots given target J and degree
#'
#' Computes the interior knot vector for a B-spline basis with a desired
#' number of basis functions J and spline degree p = degree, over the
#' domain specified by Boundary.knots = c(a, b).
#'
#' Identity used:  J = (# interior knots) + degree + 1  (with intercept = TRUE)
#' so the number of interior knots is q = J - degree - 1 (must be >= 0).
#'
#' @param J Integer, desired number of basis functions (must satisfy J >= degree + 1).
#' @param degree Integer >= 0, spline degree p.
#' @param Boundary.knots Numeric length-2, c(a, b) with a < b (the spline domain).
#' @param method "equispaced" (default) places q interior knots evenly inside (a, b);
#'   "quantile" places q interior knots at the empirical quantiles of x_grid.
#' @param x_grid Numeric vector of support points used for quantile placement
#'   (required when method = "quantile"). Quantiles used are i/(q+1), i=1..q.
#'
#' @return Numeric vector of interior knots of length q (possibly length 0 if q=0).
#' @export
make_interior_knots <- function(J,
                                degree,
                                Boundary.knots,
                                method = c("equispaced", "quantile"),
                                x_grid = NULL) {
  method <- match.arg(method)
  stopifnot(length(J) == 1L, is.finite(J), J == as.integer(J))
  stopifnot(length(degree) == 1L, is.finite(degree), degree == as.integer(degree), degree >= 0L)
  stopifnot(is.numeric(Boundary.knots), length(Boundary.knots) == 2L)
  
  a <- Boundary.knots[1]; b <- Boundary.knots[2]
  if (!(is.finite(a) && is.finite(b) && a < b)) {
    stop("Boundary.knots must be finite with a < b.")
  }
  
  q <- as.integer(J - degree - 1L)
  if (q < 0L) stop("J must satisfy J >= degree + 1 (so q = J - degree - 1 >= 0).")
  if (q == 0L) return(numeric(0L))
  
  if (method == "equispaced") {
    # Place q points strictly inside (a, b) at equal spacing
    interior <- seq(a + (b - a) / (q + 1), b - (b - a) / (q + 1), length.out = q)
  } else {
    # Quantile-based placement from x_grid
    if (is.null(x_grid) || !is.numeric(x_grid) || length(x_grid) < q + 2L) {
      stop("For method='quantile', provide numeric x_grid with length >= q + 2.")
    }
    # Ensure x_grid within [a, b] and finite
    xg <- x_grid[is.finite(x_grid) & x_grid >= a & x_grid <= b]
    if (length(xg) < q + 2L) stop("x_grid must contain enough points within [a, b].")
    probs <- (1:q) / (q + 1)
    interior <- as.numeric(stats::quantile(xg, probs = probs, names = FALSE, type = 7))
    # Guard against ties at boundaries due to discrete grids
    interior[interior <= a] <- a + .Machine$double.eps^0.5
    interior[interior >= b] <- b - .Machine$double.eps^0.5
    # Enforce strict monotonicity if tiny numerical ties occur
    if (any(diff(interior) <= 0)) {
      interior <- sort(unique(interior))
      if (length(interior) != q) {
        stop("Quantile placement produced tied knots; consider jittering x_grid or using equispaced.")
      }
    }
  }
  interior
}


#' Stepwise day basis b(m) built from B-splines (splines2 backend)
#'
#' Constructs a per-day step basis \code{B_day} (J x 365) from a continuous
#' B-spline basis defined by a knot sequence and degree. Two modes:
#' \itemize{
#'   \item \code{eval_rule = "midpoint"}: b(m) := b_spl(m - 0.5) (fast)
#'   \item \code{eval_rule = "integral"}: b(m) := \int_{m-1}^{m} b_spl(x)\,dx
#'         (via \code{splines2::ibs}, numerically stable and exact for polynomials)
#' }
#'
#' This function accepts a "full" knot vector with boundary repeats (open-uniform
#' style) \emph{or} the usual \code{splines2} interface of interior knots +
#' Boundary.knots; if only a full vector is provided, it is decomposed internally.
#'
#' @param knots Numeric, nondecreasing knot vector. Either:
#'   \itemize{
#'     \item Full knot vector with boundary repeats (length K), or
#'     \item Interior knots (no repeats at the boundary).
#'   }
#'   If \code{Boundary.knots} is missing, \code{knots} is treated as a full knot vector.
#' @param degree Integer >= 0, spline degree p.
#' @param J Optional integer (# basis functions). If supplied, checked against
#'   \eqn{length(full_knots) - degree - 1}.
#' @param eval_rule Either "integral" (default) or "midpoint".
#' @param Boundary.knots Optional numeric length-2 (c(left, right)). If provided,
#'   \code{knots} is interpreted as interior knots (no boundary repeats).
#'
#' @return Numeric matrix \code{B_day} of size J x 365 whose m-th column is b(m).
#'   Columns are clipped to nonnegative (within tolerance) and normalized to sum 1.
#'
#' @importFrom splines2 bSpline ibs
#' @export
build_B_day_from_spline <- function(knots,
                                    degree,
                                    J = NULL,
                                    eval_rule = c("integral", "midpoint"),
                                    Boundary.knots = NULL) {
  eval_rule <- match.arg(eval_rule)
  stopifnot(is.numeric(knots), length(knots) >= 1L)
  stopifnot(length(degree) == 1L, degree >= 0, is.finite(degree))
  
  # Detect whether 'knots' is a full knot vector (with boundary repeats)
  # or just interior knots + Boundary.knots. If Boundary.knots not given,
  # treat 'knots' as full and decompose to (interior, Boundary).
  if (is.null(Boundary.knots)) {
    full_knots <- knots
    K <- length(full_knots)
    if (K < (2 * (degree + 1) + 1)) {
      stop("Full knot vector too short for given degree. Need >= 2*(degree+1)+1.")
    }
    # Open-uniform full knots: first and last (degree+1) are boundary repeats
    Boundary.knots <- c(full_knots[degree + 1L], full_knots[K - degree])
    if (!isTRUE(all.equal(Boundary.knots[1], min(full_knots))) ||
        !isTRUE(all.equal(Boundary.knots[2], max(full_knots)))) {
      # Still fine: we trust the positions implied by the full vector
      Boundary.knots <- c(full_knots[degree + 1L], full_knots[K - degree])
    }
    # Interior knots are those between the boundary repeats
    if (K > 2L * (degree + 1L)) {
      interior <- full_knots[(degree + 2L):(K - degree - 1L)]
    } else {
      interior <- numeric(0L)
    }
    J_imp <- K - degree - 1L
  } else {
    # 'knots' is interior knots; Boundary.knots given explicitly
    interior <- knots
    if (length(Boundary.knots) != 2L) stop("Boundary.knots must be length 2.")
    # With splines2: #basis = length(interior) + degree + 1 (when intercept=TRUE)
    J_imp <- length(interior) + degree + 1L
  }
  
  if (!is.null(J) && J != J_imp) {
    stop(sprintf("J (%d) does not match implied number of basis functions (%d).",
                 J, J_imp))
  }
  J <- J_imp
  
  # Build B_day using splines2
  if (eval_rule == "midpoint") {
    # Evaluate the B-spline basis at midpoints m - 0.5 for m=1..365
    x_mid <- seq(0.5, 364.5, by = 1.0)
    Bx <- splines2::bSpline(
      x            = x_mid,
      knots        = interior,
      degree       = degree,
      Boundary.knots = Boundary.knots,
      intercept    = TRUE
    ) # dim: 365 x J
    B_day <- t(Bx) # J x 365
  } else {
    # Use integrated B-splines (IBS). For each basis j, IBS_j(x) = \int_{left}^x b_j(t) dt.
    # Then b(m) = IBS(m) - IBS(m-1) (vectorized across columns).
    x_edges <- 0:365
    I <- splines2::ibs(
      x            = x_edges,
      knots        = interior,
      degree       = degree,
      Boundary.knots = Boundary.knots,
      intercept    = TRUE
    ) # dim: 366 x J  (cumulative integrals at day edges)
    D <- diff(I, lag = 1L, differences = 1L) # 365 x J
    B_day <- t(D) # J x 365
  }
  
  # Numerical hygiene: clip tiny negatives; enforce column sums approx 1
  B_day[B_day < 0 & B_day > -1e-12] <- 0
  if (any(B_day < -1e-8)) {
    warning("Large negative basis values detected; knots/degree may be inconsistent.")
    B_day[B_day < 0] <- 0
  }
  col_sums <- colSums(B_day)
  if (any(!is.finite(col_sums) | col_sums <= 0)) {
    bad <- which(!is.finite(col_sums) | col_sums <= 0)
    stop(sprintf("Nonpositive or nonfinite column sums for b(m) at days: %s", paste(bad, collapse = ", ")))
  }
  B_day <- sweep(B_day, 2L, col_sums, "/")
  
  rownames(B_day) <- sprintf("b%02d", seq_len(J))
  colnames(B_day) <- sprintf("day%03d", seq_len(365))
  B_day
}


#' Construct the RSV model object
#'
#' This function validates and packages the global, precomputed components of
#' the RSV model into a standardized container for use in the likelihood and
#' related computations. The resulting object is intended to be treated as
#' read-only in downstream code.
#'
#' @param B_day Matrix or array containing the basis for the infection-age
#'   density \eqn{w(a; \beta)} evaluated on the age (day) grid.
#'
#' @param S_day Matrix or array containing the basis for the detection function
#'   \eqn{c(a; \eta)} evaluated on the same age (day) grid.
#'
#' @param phi Numeric matrix of precomputed quantities (e.g., day-averaged
#'   contributions) used in the likelihood. Its dimensions should be consistent
#'   with the basis for \eqn{w(a; \beta)} and the age grid.
#'
#' @param lambda_global Numeric vector giving the calendar-time circulation
#'   function \eqn{\lambda(t)} over the full calendar window of interest.
#'
#' @return An object of class \code{"rsv_model"}.
#'
make_rsv_model <- function(B_day, S_day, phi, lambda_global) {
  
  # ---- Basic type checks ----------------------------------------------------
  
  if (!is.matrix(B_day) && !is.array(B_day)) {
    stop("`B_day` must be a matrix or array.")
  }
  if (!is.numeric(B_day)) {
    stop("`B_day` must be numeric.")
  }
  
  if (!is.matrix(S_day) && !is.array(S_day)) {
    stop("`S_day` must be a matrix or array.")
  }
  if (!is.numeric(S_day)) {
    stop("`S_day` must be numeric.")
  }
  
  if (!is.matrix(phi)) {
    stop("`phi` must be a matrix.")
  }
  if (!is.numeric(phi)) {
    stop("`phi` must be numeric.")
  }
  
  if (!is.numeric(lambda_global) || is.matrix(lambda_global) || is.array(lambda_global)) {
    stop("`lambda_global` must be a numeric vector.")
  }
  
  # ---- Non-trivial sizes ----------------------------------------------------
  
  if (nrow(B_day) == 0L || ncol(B_day) == 0L) {
    stop("`B_day` must have positive numbers of rows and columns.")
  }
  if (nrow(S_day) == 0L || ncol(S_day) == 0L) {
    stop("`S_day` must have positive numbers of rows and columns.")
  }
  if (nrow(phi) == 0L || ncol(phi) == 0L) {
    stop("`phi` must have positive numbers of rows and columns.")
  }
  
  # ---- Dimension consistency checks ----------------------------------------
  
  # Same age grid (columns) for B_day, S_day, and phi
  if (ncol(B_day) != ncol(S_day)) {
    stop("`B_day` and `S_day` must have the same number of columns (same age grid).")
  }
  if (ncol(B_day) != ncol(phi)) {
    stop("`B_day` and `phi` must have the same number of columns (same age grid).")
  }
  
  # Same number of basis functions for w in B_day and phi
  if (nrow(phi) != nrow(B_day)) {
    stop("`phi` must have the same number of rows as `B_day` (same J basis functions).")
  }
  
  # ---- Sanity checks on numeric contents -----------------------------------
  
  if (any(!is.finite(lambda_global))) {
    stop("`lambda_global` contains non-finite values (NA, NaN, or Inf).")
  }
  
  if (any(!is.finite(phi))) {
    stop("`phi` contains non-finite values (NA, NaN, or Inf).")
  }
  
  # ---- Create the rsv_model object -----------------------------------------
  
  out <- list(
    B_day         = B_day,
    S_day         = S_day,
    phi           = phi,
    lambda_global = lambda_global
  )
  class(out) <- "rsv_model"
  
  return(out)
}


#' Construct second-difference penalty matrix
#'
#' Creates the discrete second-difference matrix D of size (p-2) x p.
#' For a parameter vector theta of length p, D %*% theta yields
#' the vector of second differences:
#'
#'   theta_{j+2} - 2*theta_{j+1} + theta_j
#'
#' for j = 1, ..., p-2.
#'
#' This matrix is typically used to define a quadratic roughness penalty
#'
#'   (alpha / 2) * || D theta ||^2
#'
#' @param p Integer >= 3. Length of the parameter vector.
#'
#' @return A dense numeric matrix of dimension (p-2) x p.
#'
#' @examples
#' D <- make_D2(6)
#' theta <- 1:6
#' D %*% theta  # should be zero (linear function)
#'
#' @export
make_D2 <- function(p) {
  # ---- checks ----
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p)) {
    stop("p must be a single finite numeric value.")
  }
  
  p <- as.integer(p)
  
  if (p < 3L) {
    stop("Second-difference penalty requires p >= 3.")
  }
  
  # ---- allocate matrix ----
  D <- matrix(0.0, nrow = p - 2L, ncol = p)
  
  # ---- fill 1, -2, 1 pattern ----
  for (i in seq_len(p - 2L)) {
    D[i, i]     <-  1.0
    D[i, i + 1] <- -2.0
    D[i, i + 2] <-  1.0
  }
  
  D
}


fit_true_curves <- function(model,
                            w_target,
                            c_target,
                            days,
                            normalize = TRUE,
                            opt_maxit = 2000L) {
  stopifnot(inherits(model, "rsv_model"))
  stopifnot(is.numeric(days), length(days) == 1L, days > 0)
  
  B_day <- model$B_day
  S_day <- model$S_day
  
  # ---- Dimension checks ----------------------------------------------------
  if (ncol(B_day) != days) {
    stop(sprintf(
      "fit_true_curves: ncol(B_day) = %d but days = %d.",
      ncol(B_day), days
    ))
  }
  if (ncol(S_day) != days) {
    stop(sprintf(
      "fit_true_curves: ncol(S_day) = %d but days = %d.",
      ncol(S_day), days
    ))
  }
  if (length(w_target) != days) {
    stop(sprintf(
      "fit_true_curves: w_target must have length = %d (got %d).",
      days, length(w_target)
    ))
  }
  if (length(c_target) != days) {
    stop(sprintf(
      "fit_true_curves: c_target must have length = %d (got %d).",
      days, length(c_target)
    ))
  }
  
  if (any(w_target < 0) || any(c_target < 0)) {
    stop("fit_true_curves: w_target and c_target must be nonnegative.")
  }
  
  # ---- Optional normalization of targets ----------------------------------
  if (normalize) {
    sw <- sum(w_target)
    sc <- sum(c_target)
    if (sw <= 0 || sc <= 0) {
      stop("fit_true_curves: cannot normalize targets with zero total mass.")
    }
    w_target <- w_target / sw
    c_target <- c_target / sc
  }
  
  # ---- Setup least-squares problems ---------------------------------------
  X_w <- t(B_day)  # days x J
  X_c <- t(S_day)  # days x K
  
  obj_w <- function(beta) {
    r <- as.numeric(X_w %*% beta - w_target)
    sum(r * r)
  }
  obj_c <- function(eta) {
    r <- as.numeric(X_c %*% eta - c_target)
    sum(r * r)
  }
  
  # ---- Initial values (projected LS) --------------------------------------
  beta_init <- pmax(as.numeric(qr.solve(X_w, w_target)), 0)
  eta_init  <- pmax(as.numeric(qr.solve(X_c, c_target)), 0)
  
  if (sum(beta_init) == 0) beta_init <- rep(1e-8, nrow(B_day))
  if (sum(eta_init)  == 0) eta_init  <- rep(1e-8, nrow(S_day))
  
  # ---- Optimization -------------------------------------------------------
  opt_w <- optim(
    par     = beta_init,
    fn      = obj_w,
    method  = "L-BFGS-B",
    lower   = rep(0, length(beta_init)),
    control = list(maxit = opt_maxit)
  )
  opt_c <- optim(
    par     = eta_init,
    fn      = obj_c,
    method  = "L-BFGS-B",
    lower   = rep(0, length(eta_init)),
    control = list(maxit = opt_maxit)
  )
  
  beta0 <- as.numeric(opt_w$par)
  eta0  <- as.numeric(opt_c$par)
  
  # ---- Reconstruct curves -------------------------------------------------
  w_basis <- as.numeric(crossprod(beta0, B_day))
  c_basis <- as.numeric(crossprod(eta0,  S_day))
  
  if (normalize) {
    sw <- sum(w_basis)
    sc <- sum(c_basis)
    if (sw <= 0 || sc <= 0) {
      stop("fit_true_curves: reconstructed curves have zero total mass.")
    }
    w_true <- w_basis / sw
    c_true <- c_basis / sc
    beta_true <- beta0 / sw
    eta_true  <- eta0  / sc
  } else {
    w_true <- w_basis
    c_true <- c_basis
    beta_true <- beta0
    eta_true  <- eta0
  }
  
  list(
    beta_true = beta_true,
    eta_true  = eta_true,
    w_true    = w_true,
    c_true    = c_true,
    w_target  = w_target,
    c_target  = c_target,
    optim_w   = opt_w,
    optim_c   = opt_c,
    normalize = normalize,
    days      = days
  )
}