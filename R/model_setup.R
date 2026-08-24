# Model setup
#
# Constructs numerical controls, spline bases, penalty matrices,
# model objects, and simulation truth used throughout the RSV workflow.


# Numerical controls -------------------------------------------------------


#' Construct numerical controls for RSV computations
#'
#' Creates the numerical settings used to validate and stabilize intermediate
#' quantities throughout the RSV probability, likelihood, and estimation
#' workflow. These settings control numerical behavior without changing the
#' statistical model.
#'
#' @param check_bounds Logical scalar indicating whether functions should stop
#'   when intermediate quantities exceed their theoretical bounds beyond the
#'   permitted numerical tolerance.
#' @param warn_on_clip Logical scalar indicating whether to issue warnings when
#'   numerical clipping or stabilization is applied.
#' @param eps_tau Positive numeric scalar used as the lower bound for no-visit
#'   probabilities and related quantities that may appear inside logarithms.
#' @param tol_clip Nonnegative numeric scalar used to distinguish small
#'   floating-point deviations from substantive bound violations.
#' @param eps_log Positive numeric scalar used as the lower bound for
#'   nonnegative quantities before applying logarithms.
#'
#' @return An `rsv_control` object containing:
#' \describe{
#'   \item{\code{check_bounds}}{Logical flag for strict bound checking.}
#'   \item{\code{warn_on_clip}}{Logical flag for warnings when clipping or
#'     stabilization occurs.}
#'   \item{\code{eps_tau}}{Positive lower bound for no-visit probabilities and
#'     related quantities.}
#'   \item{\code{tol_clip}}{Nonnegative tolerance for numerical bound
#'     violations.}
#'   \item{\code{eps_log}}{Positive lower bound for quantities evaluated inside
#'     logarithms.}
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


# Simulation model construction -------------------------------------------


#' Construct an RSV simulation model from a circulation curve
#'
#' Builds the spline bases, penalty matrices, and `rsv_model` object used in
#' simulations with a user-supplied seasonal RSV circulation curve. The
#' circulation curve is defined on the calendar-day grid, while the spline
#' bases represent the infection-age and healthcare-visit curves on the daily
#' age grid.
#'
#' @param lambda_global Nonnegative numeric vector of length `calendar_len`
#'   containing the seasonal RSV circulation curve \eqn{\lambda(t)} on the
#'   calendar-day grid. Values must be finite.
#' @param J Positive integer scalar giving the number of basis functions for
#'   the infection-age curve and the length of \eqn{\beta}.
#' @param K Positive integer scalar giving the number of basis functions for
#'   the healthcare-visit curve and the length of \eqn{\eta}.
#' @param degree Nonnegative integer scalar giving the degree of the B-spline
#'   bases.
#' @param days Positive integer scalar giving the number of days in the age
#'   grid used to construct `B_day`, `S_day`, and `phi`. In the current
#'   implementation, this must be 365 and equal `age_len`.
#' @param age_len Positive integer scalar giving the length of the modeled age
#'   window. In the current implementation, this must be 365 and equal `days`.
#' @param birth_max Positive integer scalar giving the largest one-based birth
#'   index used in the simulation.
#' @param calendar_len Positive integer scalar giving the number of days in the
#'   calendar grid. It must equal `length(lambda_global)` and be at least
#'   `birth_max + age_len`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{model}}{An `rsv_model` object containing the daily spline
#'     bases, seasonal circulation curve, second-difference matrices, and
#'     penalty matrices.}
#'   \item{\code{B_day}}{Numeric `J x 365` matrix containing the daily basis
#'     for the infection-age curve.}
#'   \item{\code{S_day}}{Numeric `K x 365` matrix containing the daily basis
#'     for the healthcare-visit curve.}
#'   \item{\code{phi}}{Numeric `J x 365` matrix containing the day-integrated
#'     infection-age basis used to construct subject-level cumulative kernels.}
#'   \item{\code{lambda_global}}{Nonnegative numeric vector containing the
#'     supplied seasonal RSV circulation curve.}
#'   \item{\code{calendar_len}}{Integer scalar giving the calendar-grid
#'     length.}
#'   \item{\code{birth_max}}{Integer scalar giving the largest permitted birth
#'     index.}
#'   \item{\code{age_len}}{Integer scalar giving the modeled age-window
#'     length.}
#'   \item{\code{lambda_params}}{Named list identifying the circulation curve
#'     as a user-provided custom input.}
#' }
#'
#' @details
#' For a subject with birth index \eqn{B_i}, age day \eqn{a} corresponds to
#' calendar index \eqn{B_i + a}. The calendar grid must therefore extend
#' through at least `birth_max + age_len`.
#'
#' The model includes second-difference penalty matrices
#' \eqn{M_\beta = D_\beta^\top D_\beta} and
#' \eqn{M_\eta = D_\eta^\top D_\eta} for the infection-age and
#' healthcare-visit coefficient vectors, respectively.
make_sim_model_lambda_custom <- function(lambda_global,
                                         J = 18L,
                                         K = 18L,
                                         degree = 3L,
                                         days = 365L,
                                         age_len = 365L,
                                         birth_max = 365L,
                                         calendar_len = 730L) {
  # ---- Input validation ----
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


# Spline basis construction -----------------------------------------------


#' Construct day-integrated spline bases for RSV simulations
#'
#' Builds the spline basis matrices used to represent the infection-age curve,
#' the healthcare-visit curve, and the cumulative infection kernel on the
#' daily age grid. Each daily basis vector is obtained by integrating a
#' continuous B-spline basis over the corresponding day interval.
#'
#' @param J Positive integer scalar giving the number of basis functions for
#'   the infection-age curve and the number of rows in `B_day` and `phi`. It
#'   must be at least `degree + 1`.
#' @param K Positive integer scalar giving the number of basis functions for
#'   the healthcare-visit curve and the number of rows in `S_day`. It must be
#'   at least `degree + 1`.
#' @param degree Nonnegative integer scalar giving the degree of the B-spline
#'   bases.
#' @param days Positive integer scalar giving the number of daily intervals in
#'   the age grid. In the current implementation, this must be 365.
#' @param knot_method Character scalar specifying the interior-knot placement
#'   method. Supported values are `"quantile"` and `"equispaced"`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{B_day}}{Numeric `J x 365` matrix containing the
#'     day-integrated basis for the infection-age curve. Column \eqn{m}
#'     represents the basis integrated over \eqn{(m-1,m]}.}
#'   \item{\code{S_day}}{Numeric `K x 365` matrix containing the
#'     day-integrated basis for the healthcare-visit curve. Column \eqn{m}
#'     represents the basis integrated over \eqn{(m-1,m]}.}
#'   \item{\code{phi}}{Numeric `J x 365` matrix containing the
#'     day-integrated infection-age basis used in cumulative kernel
#'     calculations.}
#'   \item{\code{Boundary.knots}}{Numeric vector of length 2 containing the
#'     spline-domain endpoints `c(0, days)`.}
#'   \item{\code{degree}}{Nonnegative integer scalar giving the spline degree.}
#'   \item{\code{days}}{Positive integer scalar giving the number of daily
#'     intervals.}
#' }
#'
#' @details
#' For day \eqn{m}, the infection-age basis is
#' \deqn{B_{\mathrm{day}}[,m] = \int_{m-1}^{m} b(x)\,dx,}
#' with an analogous construction for the healthcare-visit basis.
#'
#' When `K == J`, `S_day` is set equal to `B_day`. Otherwise, a separate
#' `K`-dimensional basis is constructed. In the current implementation,
#' `B_day` and `phi` contain the same integrated infection-age basis but are
#' retained as separate components because they serve different downstream
#' roles.
make_sim_bases <- function(J = 18L,
                           K = 18L,
                           degree = 3L,
                           days = 365L,
                           knot_method = c("quantile", "equispaced")) {
  knot_method <- match.arg(knot_method)
  stopifnot(J > 0L, K > 0L, days > 0L)
  
  day_grid <- seq_len(days)
  Boundary.knots <- c(0, days)  # domain for the spline basis
  
  # ---- Infection-age basis -----------------------------------------------
  knots_w <- make_interior_knots(
    J              = J,
    degree         = degree,
    Boundary.knots = Boundary.knots,
    method         = knot_method,
    x_grid         = day_grid
  )
  
  B_day <- build_B_day_from_spline(
    knots          = knots_w,
    degree         = degree,
    J              = J,
    eval_rule      = "integral",
    Boundary.knots = Boundary.knots
  )
  
  # ---- Cumulative-kernel basis -------------------------------------------
  phi <- build_B_day_from_spline(
    knots          = knots_w,
    degree         = degree,
    J              = J,
    eval_rule      = "integral",   
    Boundary.knots = Boundary.knots
  )
  
  # ---- Healthcare-visit basis --------------------------------------------
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
  
  # ---- Dimension checks --------------------------------------------------
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


#' Construct interior knots for a B-spline basis
#'
#' Determines the interior knots required to construct a B-spline basis with
#' `J` basis functions and the specified spline degree. Knots may be placed at
#' equally spaced locations or at empirical quantiles of a supplied grid.
#'
#' @param J Positive integer scalar giving the desired number of basis
#'   functions. It must be at least `degree + 1`.
#' @param degree Nonnegative integer scalar giving the spline degree.
#' @param Boundary.knots Numeric vector of length 2 containing the lower and
#'   upper endpoints of the spline domain. Both values must be finite, and the
#'   lower endpoint must be smaller than the upper endpoint.
#' @param method Character scalar specifying the knot-placement method.
#'   Supported values are `"equispaced"` and `"quantile"`.
#' @param x_grid Numeric vector containing the support points used for
#'   quantile-based placement. It is required when `method = "quantile"` and
#'   must contain enough finite values within `Boundary.knots`.
#'
#' @return Numeric vector of length
#'   \eqn{q = J - \mathrm{degree} - 1} containing the interior knots. The
#'   returned vector has length zero when \eqn{q = 0}.
#'
#' @details
#' With an intercept included, the number of B-spline basis functions satisfies
#' \deqn{J = q + \mathrm{degree} + 1,}
#' where \eqn{q} is the number of interior knots.
#'
#' For equispaced placement, the knots divide the interior of the spline domain
#' into equal intervals. For quantile-based placement, the knots are placed at
#' probabilities \eqn{r/(q+1)} for \eqn{r = 1,\ldots,q} using the finite values
#' of `x_grid` that fall within the spline domain.
#'
#' Quantile-based knots are adjusted away from the boundary when necessary.
#' The function stops if quantile placement produces tied interior knots.
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
    if (is.null(x_grid) || !is.numeric(x_grid) || length(x_grid) < q + 2L) {
      stop("For method='quantile', provide numeric x_grid with length >= q + 2.")
    }
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


#' Construct a daily B-spline basis matrix
#'
#' Constructs a B-spline basis on the 365-day age grid using either basis
#' values at daily midpoints or basis integrals over daily intervals. The
#' resulting columns are normalized to sum to one.
#'
#' @param knots Numeric vector specifying the spline knots. When
#'   `Boundary.knots` is `NULL`, this must be a full knot vector containing
#'   repeated boundary knots. When `Boundary.knots` is supplied, this vector is
#'   interpreted as the interior knots.
#' @param degree Nonnegative integer scalar giving the spline degree.
#' @param J `NULL` or a positive integer scalar giving the expected number of
#'   basis functions. When supplied, it must match the number implied by the
#'   knots and degree.
#' @param eval_rule Character scalar specifying how to construct each daily
#'   basis vector. Supported values are `"integral"` and `"midpoint"`.
#' @param Boundary.knots `NULL` or a numeric vector of length 2 containing the
#'   lower and upper endpoints of the spline domain. When supplied, `knots` is
#'   interpreted as a vector of interior knots.
#'
#' @return Numeric `J x 365` matrix containing the daily B-spline basis.
#'   Column \eqn{m} corresponds to age day \eqn{m}, and each column is
#'   normalized to sum to one.
#'
#' @details
#' When `eval_rule = "integral"`, column \eqn{m} is constructed as
#' \deqn{B_{\mathrm{day}}[,m]
#'   = \int_{m-1}^{m} b(x)\,dx.}
#'
#' When `eval_rule = "midpoint"`, column \eqn{m} is constructed by evaluating
#' the basis at the midpoint \eqn{m - 1/2}.
#'
#' If a full knot vector is supplied, the first and last `degree + 1` entries
#' are treated as repeated boundary knots. Otherwise, the number of basis
#' functions is
#' \deqn{J = q + \mathrm{degree} + 1,}
#' where \eqn{q} is the number of interior knots.
#'
#' Small negative values caused by floating-point error are set to zero.
#' Larger negative values produce a warning and are also set to zero before
#' column normalization.
build_B_day_from_spline <- function(knots,
                                    degree,
                                    J = NULL,
                                    eval_rule = c("integral", "midpoint"),
                                    Boundary.knots = NULL) {
  eval_rule <- match.arg(eval_rule)
  stopifnot(is.numeric(knots), length(knots) >= 1L)
  stopifnot(length(degree) == 1L, degree >= 0, is.finite(degree))
  
  # ---- Resolve knot representation -----------------------------------------
  if (is.null(Boundary.knots)) {
    full_knots <- knots
    K <- length(full_knots)
    if (K < (2 * (degree + 1) + 1)) {
      stop("Full knot vector too short for given degree. Need >= 2*(degree+1)+1.")
    }
    # The first and last degree + 1 knots define the repeated boundaries.
    Boundary.knots <- c(full_knots[degree + 1L], full_knots[K - degree])
    if (!isTRUE(all.equal(Boundary.knots[1], min(full_knots))) ||
        !isTRUE(all.equal(Boundary.knots[2], max(full_knots)))) {
      Boundary.knots <- c(full_knots[degree + 1L], full_knots[K - degree])
    }
    if (K > 2L * (degree + 1L)) {
      interior <- full_knots[(degree + 2L):(K - degree - 1L)]
    } else {
      interior <- numeric(0L)
    }
    J_imp <- K - degree - 1L
  } else {
    interior <- knots
    if (length(Boundary.knots) != 2L) stop("Boundary.knots must be length 2.")
    # With an intercept, J = number of interior knots + degree + 1.
    J_imp <- length(interior) + degree + 1L
  }
  
  if (!is.null(J) && J != J_imp) {
    stop(sprintf("J (%d) does not match implied number of basis functions (%d).",
                 J, J_imp))
  }
  J <- J_imp
  
  # ---- Construct daily basis -----------------------------------------------
  if (eval_rule == "midpoint") {
    # Evaluate the basis at the midpoint of each daily interval.
    x_mid <- seq(0.5, 364.5, by = 1.0)
    Bx <- splines2::bSpline(
      x            = x_mid,
      knots        = interior,
      degree       = degree,
      Boundary.knots = Boundary.knots,
      intercept    = TRUE
    ) 
    B_day <- t(Bx) 
  } else {
    # Obtain each daily basis vector by differentiating cumulative integrals
    # evaluated at consecutive day boundaries.
    x_edges <- 0:365
    I <- splines2::ibs(
      x            = x_edges,
      knots        = interior,
      degree       = degree,
      Boundary.knots = Boundary.knots,
      intercept    = TRUE
    ) 
    D <- diff(I, lag = 1L, differences = 1L) 
    B_day <- t(D) 
  }
  
  # ---- Stabilize and normalize ---------------------------------------------
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


#' Construct an RSV model object
#'
#' Validates and combines the spline basis matrices, cumulative-kernel basis,
#' and calendar-time RSV circulation curve into the model object used
#' throughout the probability, likelihood, estimation, and simulation
#' workflow.
#'
#' @param B_day Numeric `J x D` matrix or array containing the daily basis for
#'   the infection-age curve `w(a)`.
#' @param S_day Numeric `K x D` matrix or array containing the daily basis for
#'   the healthcare-visit curve `c(a)`.
#' @param phi Numeric `J x D` matrix containing the day-integrated
#'   infection-age basis used to construct subject-specific cumulative kernels.
#' @param lambda_global Numeric vector containing the RSV circulation curve
#'   \eqn{\lambda(t)} on the full calendar-time grid.
#'
#' @return An `rsv_model` object containing:
#' \describe{
#'   \item{\code{B_day}}{Numeric `J x D` matrix or array containing the
#'     infection-age basis.}
#'   \item{\code{S_day}}{Numeric `K x D` matrix or array containing the
#'     healthcare-visit basis.}
#'   \item{\code{phi}}{Numeric `J x D` matrix containing the day-integrated
#'     infection-age basis used in cumulative-kernel calculations.}
#'   \item{\code{lambda_global}}{Numeric vector containing the calendar-time
#'     RSV circulation curve.}
#' }
#'
#' @details
#' `B_day`, `S_day`, and `phi` must have the same number of columns so that
#' they represent the same age grid. `phi` must also have the same number of
#' rows as `B_day` because both use the same `J` infection-age basis
#' functions.
make_rsv_model <- function(B_day, S_day, phi, lambda_global) {
  
  # ---- Input validation ----------------------------------------------------
  
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
  
  # ---- Validate dimensions -------------------------------------------------
  
  if (nrow(B_day) == 0L || ncol(B_day) == 0L) {
    stop("`B_day` must have positive numbers of rows and columns.")
  }
  if (nrow(S_day) == 0L || ncol(S_day) == 0L) {
    stop("`S_day` must have positive numbers of rows and columns.")
  }
  if (nrow(phi) == 0L || ncol(phi) == 0L) {
    stop("`phi` must have positive numbers of rows and columns.")
  }
  
  # ---- Dimension consistency -----------------------------------------------
  
  if (ncol(B_day) != ncol(S_day)) {
    stop("`B_day` and `S_day` must have the same number of columns (same age grid).")
  }
  if (ncol(B_day) != ncol(phi)) {
    stop("`B_day` and `phi` must have the same number of columns (same age grid).")
  }
  
  if (nrow(phi) != nrow(B_day)) {
    stop("`phi` must have the same number of rows as `B_day` (same J basis functions).")
  }
  
  # ---- Validate numeric contents -------------------------------------------
  
  if (any(!is.finite(lambda_global))) {
    stop("`lambda_global` contains non-finite values (NA, NaN, or Inf).")
  }
  
  if (any(!is.finite(phi))) {
    stop("`phi` contains non-finite values (NA, NaN, or Inf).")
  }
  
  # ---- Construct model object ----------------------------------------------
  
  out <- list(
    B_day         = B_day,
    S_day         = S_day,
    phi           = phi,
    lambda_global = lambda_global
  )
  class(out) <- "rsv_model"
  
  return(out)
}


#' Construct a second-difference matrix
#'
#' Constructs the discrete second-difference matrix used to penalize roughness
#' in spline coefficient vectors. For a coefficient vector \eqn{\theta} of
#' length \eqn{p}, the product \eqn{D\theta} contains the second differences
#'
#' \deqn{
#'   \theta_{j+2} - 2\theta_{j+1} + \theta_j,
#'   \qquad j = 1,\ldots,p-2.
#' }
#'
#' @param p Integer scalar giving the length of the coefficient vector. It must
#'   be at least 3.
#'
#' @return Numeric `(p - 2) x p` matrix whose rows contain the
#'   second-difference pattern \eqn{(1,-2,1)}.
#'
#' @details
#' If \eqn{M = D^\top D}, then the quadratic roughness penalty can be written as
#'
#' \deqn{
#'   \frac{\alpha}{2}\theta^\top M\theta
#'   =
#'   \frac{\alpha}{2}\lVert D\theta\rVert^2.
#' }
#'
#' This penalty discourages large changes in successive coefficient slopes and
#' therefore favors smoother spline curves.
make_D2 <- function(p) {
  # ---- Input validation ----------------------------------------------------
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p)) {
    stop("p must be a single finite numeric value.")
  }
  
  p <- as.integer(p)
  
  if (p < 3L) {
    stop("Second-difference penalty requires p >= 3.")
  }
  
  D <- matrix(0.0, nrow = p - 2L, ncol = p)
  
  # Fill each row with the second-difference pattern (1, -2, 1).
  for (i in seq_len(p - 2L)) {
    D[i, i]     <-  1.0
    D[i, i + 1] <- -2.0
    D[i, i + 2] <-  1.0
  }
  
  D
}


# Simulation truth construction -------------------------------------------


#' Project target curves onto the model spline bases
#'
#' Approximates target infection-age and healthcare-visit curves using the
#' spline bases stored in an `rsv_model` object. The resulting nonnegative
#' coefficients and reconstructed curves define the simulation truth used in
#' the reference workflow.
#'
#' @param model An `rsv_model` object containing `B_day` and `S_day`.
#' @param w_target Nonnegative numeric vector of length `days` containing the
#'   target infection-age curve.
#' @param c_target Nonnegative numeric vector of length `days` containing the
#'   target healthcare-visit curve.
#' @param days Positive integer scalar giving the age-grid length.
#' @param normalize Logical scalar indicating whether the target and
#'   reconstructed curves should be normalized to sum to one.
#' @param opt_maxit Positive integer scalar giving the maximum number of
#'   optimization iterations.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{beta_true}}{Nonnegative infection-age spline coefficients.}
#'   \item{\code{eta_true}}{Nonnegative healthcare-visit spline coefficients.}
#'   \item{\code{w_true}}{Model-represented infection-age curve.}
#'   \item{\code{c_true}}{Model-represented healthcare-visit curve.}
#'   \item{\code{w_target}}{Target infection-age curve used for projection.}
#'   \item{\code{c_target}}{Target healthcare-visit curve used for projection.}
#'   \item{\code{optim_w}}{Optimization result for the infection-age curve.}
#'   \item{\code{optim_c}}{Optimization result for the healthcare-visit curve.}
#'   \item{\code{normalize}}{Logical scalar indicating whether normalization
#'     was applied.}
#'   \item{\code{days}}{Integer scalar giving the age-grid length.}
#' }
#'
#' @details
#' Each target is approximated by constrained least squares with nonnegative
#' spline coefficients. The returned `w_true` and `c_true` are the resulting
#' model-represented curves and may differ slightly from the supplied targets.
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
  
  # ---- Normalize targets ---------------------------------------------------
  if (normalize) {
    sw <- sum(w_target)
    sc <- sum(c_target)
    if (sw <= 0 || sc <= 0) {
      stop("fit_true_curves: cannot normalize targets with zero total mass.")
    }
    w_target <- w_target / sw
    c_target <- c_target / sc
  }
  
  # ---- Setup least-squares projections ------------------------------------
  X_w <- t(B_day)
  X_c <- t(S_day)
  
  obj_w <- function(beta) {
    r <- as.numeric(X_w %*% beta - w_target)
    sum(r * r)
  }
  obj_c <- function(eta) {
    r <- as.numeric(X_c %*% eta - c_target)
    sum(r * r)
  }
  
  # ---- Construct initial values --------------------------------------------
  beta_init <- pmax(as.numeric(qr.solve(X_w, w_target)), 0)
  eta_init  <- pmax(as.numeric(qr.solve(X_c, c_target)), 0)
  
  if (sum(beta_init) == 0) beta_init <- rep(1e-8, nrow(B_day))
  if (sum(eta_init)  == 0) eta_init  <- rep(1e-8, nrow(S_day))
  
  # ---- Constrained optimization --------------------------------------------
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
  
  # ---- Reconstruct model curves --------------------------------------------
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