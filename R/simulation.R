# Simulation
#
# Generates synthetic subjects, constructs subject-specific visit
# distributions, and samples healthcare-visit outcomes for the RSV workflow.


# Subject generation -------------------------------------------------------


#' Generate subjects for RSV simulation
#'
#' Constructs a synthetic cohort before healthcare visits are simulated. Each
#' subject is assigned a birth index sampled uniformly from the permitted
#' birth-index range and an initial missing healthcare-visit age.
#'
#' @param n Positive integer scalar giving the number of subjects to generate.
#' @param model An `rsv_model` object used to check simulation dimensions when
#'   an `age_len` component is available.
#' @param birth_max Positive integer scalar giving the largest birth index that
#'   can be sampled. This argument must be supplied. If `model$age_len` exists,
#'   `birth_max` must equal `model$age_len`.
#' @param seed Integer scalar or `NULL`. If supplied, it is used to make birth
#'   sampling reproducible, and the previous random-number-generator state is
#'   restored on exit.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{data}}{An `rsv_data` object containing the generated subjects.}
#'   \item{\code{subjects}}{List of `n` subject records created by
#'     `make_subject()`.}
#'   \item{\code{birth_days}}{Integer vector of length `n` containing the
#'     sampled birth indices.}
#'   \item{\code{birth_max}}{Positive integer scalar giving the largest birth
#'     index used for sampling.}
#' }
#'
#' @details
#' Birth indices are sampled independently and uniformly with replacement from
#' `1:birth_max`. Subject identifiers are assigned sequentially from 1 through
#' `n`, `visit_age` is initialized to `NA_integer_`, and the default likelihood
#' weight from `make_subject()` is one.
make_sim_subjects <- function(n, model, birth_max = NULL, seed = NULL) {
  
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1)
  
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed)
  }
  
  # Determine the permitted birth-index range.
  has_age_len <- !is.null(model$age_len)
  
  if (!has_age_len) {
    # Require birth_max when the model does not define age_len.
    if (is.null(birth_max)) {
      stop("make_sim_subjects: model has no 'age_len', so you must provide 'birth_max'.")
    }
  } else {
    # Require birth_max to match the model age window when age_len is defined.
    if (is.null(birth_max)) {
      stop("make_sim_subjects: model has 'age_len', so you must provide 'birth_max' and it must equal model$age_len.")
    }
    if (as.integer(model$age_len) != as.integer(birth_max)) {
      stop(sprintf(
        "make_sim_subjects: mismatch: model$age_len=%d but birth_max=%d. They must be equal.",
        as.integer(model$age_len), as.integer(birth_max)
      ))
    }
  }
  
  birth_max <- as.integer(birth_max)
  if (birth_max <= 0L) stop("make_sim_subjects: 'birth_max' must be a positive integer.")
  
  # Sample birth indices and construct subject records.
  birth_days <- sample.int(
    birth_max,
    size = as.integer(n),
    replace = TRUE
    )
  
  subject_list <- vector("list", length = as.integer(n))
  for (i in seq_len(as.integer(n))) {
    subject_list[[i]] <- make_subject(
      id         = i,
      birth_index = birth_days[i],
      visit_age  = NA_integer_
    )
  }
  
  data <- make_rsv_data(subjects = subject_list)
  
  list(
    data       = data,
    subjects   = subject_list,
    birth_days = birth_days,
    birth_max  = birth_max
  )
}


#' Construct an RSV subject record
#'
#' Creates a standardized subject record containing the birth timing,
#' healthcare-visit outcome, and likelihood weight used throughout the RSV
#' workflow.
#'
#' @param id Optional subject identifier.
#' @param birth_index Integer scalar giving the subject's one-based calendar
#'   birth index. This index is used to align the seasonal RSV circulation
#'   curve with the subject's age.
#' @param visit_age Integer scalar in 1:365 giving the age in days of the
#'   healthcare visit, or `NA` if no visit occurred during the first year.
#' @param weight Numeric scalar giving the multiplicative weight applied to
#'   the subject's likelihood contribution. The default is 1.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{id}}{Optional subject identifier.}
#'   \item{\code{birth_index}}{Integer scalar giving the one-based calendar
#'     birth index.}
#'   \item{\code{visit_age}}{Integer scalar in 1:365 giving the
#'     healthcare-visit age, or `NA_integer_` if no visit occurred during the
#'     first year.}
#'   \item{\code{weight}}{Numeric scalar giving the subject's likelihood
#'     weight.}
#' }
make_subject <- function(id = NULL, birth_index, visit_age, weight = 1) {
  
  if (length(birth_index) != 1L || !is.finite(birth_index)) {
    stop("birth_index must be a single finite number.")
  }
  
  # Validate the visit age only when a visit is observed.
  if (!is.na(visit_age)) {
    if (visit_age < 1 || visit_age > 365) {
      stop("visit_age must be in 1:365 or NA.")
    }
  }
  
  
  list(
    id = id,
    birth_index = as.integer(birth_index),
    visit_age = if (is.na(visit_age)) NA_integer_ else as.integer(visit_age),
    weight = weight
  )
}


#' Construct an RSV data object
#'
#' Validates and packages a list of subject records into the standardized data
#' container used throughout the RSV workflow.
#'
#' @param subjects List of subject records. Each element must be a list
#'   containing at least `birth_index`, `visit_age`, and `weight`. Additional
#'   fields, such as `id`, are allowed.
#'
#' @return An `rsv_data` object containing:
#' \describe{
#'   \item{\code{subjects}}{The validated list of subject records.}
#' }
make_rsv_data <- function(subjects) {
  
  if (!is.list(subjects)) {
    stop("`subjects` must be a list of subject objects.")
  }
  
  # Allow an empty subject list, but warn because it is atypical.
  if (length(subjects) == 0L) {
    warning("`subjects` list is empty. Likelihood will be zero.")
  }
  
  required_fields <- c("birth_index", "visit_age", "weight")
  
  for (i in seq_along(subjects)) {
    
    subj <- subjects[[i]]
    
    if (!is.list(subj)) {
      stop(sprintf("subjects[[%d]] is not a list.", i))
    }
    
    missing_fields <- setdiff(required_fields, names(subj))
    if (length(missing_fields) > 0L) {
      stop(sprintf(
        "subjects[[%d]] is missing required fields: %s",
        i,
        paste(missing_fields, collapse = ", ")
      ))
    }
    
    # Apply lightweight type checks to required numeric fields.
    if (!is.null(subj$birth_index) && !is.numeric(subj$birth_index)) {
      stop(sprintf("subjects[[%d]]$birth_index must be numeric/integer.", i))
    }
    if (!is.null(subj$weight) && !is.numeric(subj$weight)) {
      stop(sprintf("subjects[[%d]]$weight must be numeric.", i))
    }
  }
  
  out <- list(subjects = subjects)
  class(out) <- "rsv_data"
  
  return(out)
}


# Visit distribution construction -----------------------------------------


#' Construct visit distributions for all subjects
#'
#' Constructs one subject-specific distribution over healthcare-visit age and
#' the no-visit outcome for each set of subject-level precomputations.
#'
#' @param subj_pre_list Non-empty list of subject-level precomputations. Each
#'   element must contain `lambda_i`, a numeric vector of length 365, and
#'   `V_i`, a numeric J x 366 matrix.
#' @param w_vec Numeric vector of length 365 containing the infection-age
#'   curve over the daily age grid.
#' @param c_vec Numeric vector of length 365 containing the healthcare-visit
#'   curve over the daily age grid.
#' @param beta Numeric vector of length J containing the infection-age spline
#'   coefficients.
#' @param eta Numeric vector of length K containing the healthcare-visit
#'   spline coefficients. This argument is passed through to `pmf_i()`.
#' @param S_day Numeric K x 365 matrix containing the healthcare-visit spline
#'   basis. This argument is passed through to `pmf_i()`.
#' @param control An `rsv_control` object controlling numerical checks and
#'   validation behavior.
#' @param tol_prob Positive numeric scalar giving the tolerance used when
#'   checking that each subject distribution sums to one.
#'
#' @return A list with one named numeric vector of length 366 per subject.
#'   Each vector contains the masses for visit days 1 through 365 followed by
#'   the no-visit mass.
#'
#' @details
#' Each element of `subj_pre_list` is passed to `pmf_i()` together with the
#' supplied curve values, coefficients, basis, and numerical controls. The
#' returned vectors are named `day001` through `day365` and `no_visit`.
pmf_dataset <- function(subj_pre_list,
                        w_vec,
                        c_vec,
                        beta,
                        eta,
                        S_day,
                        control,
                        tol_prob = 1e-8) {
  
  if (!is.list(subj_pre_list) || length(subj_pre_list) == 0L) {
    stop("subj_pre_list must be a non-empty list.")
  }
  
  pmf_list <- lapply(subj_pre_list, function(subj_pre) {
    
    if (is.null(subj_pre$lambda_i) || is.null(subj_pre$V_i)) {
      stop("Each element of subj_pre_list must contain lambda_i and V_i.")
    }
    
    pmf_i(
      lambda_i = subj_pre$lambda_i,
      V_i      = subj_pre$V_i,
      w_vec    = w_vec,
      c_vec    = c_vec,
      beta     = beta,
      eta      = eta,
      S_day    = S_day,
      control  = control, 
      tol_prob = tol_prob
    )
  })
  
  pmf_list
}


#' Construct a visit distribution for one subject
#'
#' Constructs the subject-specific distribution over healthcare-visit age and
#' the no-visit outcome on the 365-day age grid.
#'
#' @param lambda_i Numeric vector of length 365 containing the
#'   subject-specific RSV circulation curve.
#' @param V_i Numeric J x 366 matrix containing the cumulative subject-level
#'   kernels. Column 1 represents day zero.
#' @param w_vec Numeric vector of length 365 containing the infection-age
#'   curve over the daily age grid.
#' @param c_vec Numeric vector of length 365 containing the healthcare-visit
#'   curve over the daily age grid.
#' @param beta Numeric vector of length J containing the infection-age spline
#'   coefficients.
#' @param eta Numeric vector of length K containing the healthcare-visit
#'   spline coefficients. Retained for interface consistency but not used
#'   directly in the current implementation.
#' @param S_day Numeric K x 365 matrix containing the healthcare-visit spline
#'   basis. Retained for interface consistency but not used directly in the
#'   current implementation.
#' @param control An `rsv_control` object controlling the optional
#'   normalization check.
#' @param tol_prob Positive numeric scalar giving the tolerance for the
#'   optional sum-to-one check.
#'
#' @return Named numeric vector of length 366 containing the visit-day masses
#'   for days 1 through 365 followed by the no-visit mass. The entries are
#'   named `day001` through `day365` and `no_visit`.
#'
#' @details
#' The visit density is evaluated at each of the 365 daily grid points by
#' `mass_i_visit()`, and these values are used as visit-day masses for
#' simulation. The no-visit mass is then defined as the remaining mass,
#' \deqn{
#'   p_i(\mathrm{no\ visit})
#'   =
#'   1 - \sum_{m=1}^{365} p_i(\mathrm{visit\ on\ day}\ m).
#' }
#' This complement construction makes normalization part of the simulation
#' distribution itself. When `control$check_bounds` is `TRUE`, the function
#' verifies that the resulting 366 entries sum to one within `tol_prob`.
pmf_i <- function(lambda_i,
                  V_i,
                  w_vec,
                  c_vec,
                  beta,
                  eta,
                  S_day,
                  control, 
                  tol_prob = 1e-8) {
  
  # Compute the visit-day masses on the daily grid.
  visit_mass <- mass_i_visit(
    lambda_i = lambda_i,
    V_i      = V_i,
    w_vec    = w_vec,
    c_vec    = c_vec,
    beta     = beta,
    control  = control
  )
  
  # Assign the remaining probability mass to the no-visit outcome.
  tau_mass <- 1 - sum(visit_mass)
  
  pmf <- c(visit_mass, tau_mass)
  
  names(pmf) <- c(
    sprintf("day%03d", seq_len(365)),
    "no_visit"
  )
  
  # Check normalization when bound checking is enabled.
  if (isTRUE(control$check_bounds)) {
    
    total <- sum(pmf)
    
    if (!is.finite(total) || abs(total - 1) > tol_prob) {
      stop(sprintf(
        "pmf_i does not sum to 1 (sum = %.12f).",
        total
      ))
    }
  }
  
  pmf
}


#' Compute visit-day masses
#'
#' Evaluates the visit density at the 365 daily grid points used to construct
#' a subject's simulation distribution. These values are treated as visit-day
#' masses after exponentiating the output of `logmass_i_visit()`.
#'
#' @param lambda_i Numeric vector of length 365 containing the
#'   subject-specific RSV circulation curve.
#' @param V_i Numeric J x 366 matrix containing the cumulative subject-level
#'   kernels. Column 1 represents day zero.
#' @param w_vec Numeric vector of length 365 containing the infection-age
#'   curve over the daily age grid.
#' @param c_vec Numeric vector of length 365 containing the healthcare-visit
#'   curve over the daily age grid.
#' @param beta Numeric vector of length J containing the infection-age spline
#'   coefficients.
#' @param control An `rsv_control` object retained for interface consistency
#'   but not used directly in the current implementation.
#'
#' @return Numeric vector of length 365 containing the visit-day masses.
#'
#' @details
#' For age day \eqn{m}, the visit quantity is
#' \deqn{
#'   f_i^{\mathrm{visit}}(m)
#'   =
#'   c(m)\lambda_i(m)w(m)\exp\{-H_i(m)\}.
#' }
#' These daily-grid values are used as visit-day masses by `pmf_i()`.
mass_i_visit <- function(lambda_i,
                         V_i,
                         w_vec,
                         c_vec,
                         beta,
                         control) {
  
  log_mass_vec <- logmass_i_visit(
    lambda_i = lambda_i,
    V_i      = V_i,
    w_vec    = w_vec,
    c_vec    = c_vec,
    beta     = beta,
    control  = control
  )
  
  exp(log_mass_vec)
}


#' Compute log visit-day masses
#'
#' Computes the log visit quantity at each of the 365 daily grid points used
#' to construct a subject's simulation distribution.
#'
#' @param lambda_i Numeric vector of length 365 containing the strictly
#'   positive subject-specific RSV circulation curve.
#' @param V_i Numeric J x 366 matrix containing the cumulative subject-level
#'   kernels. Column 1 represents day zero, and column `m + 1` represents
#'   \eqn{v_i(m)}.
#' @param w_vec Numeric vector of length 365 containing the strictly positive
#'   infection-age curve over the daily age grid.
#' @param c_vec Numeric vector of length 365 containing the strictly positive
#'   healthcare-visit curve over the daily age grid.
#' @param beta Numeric vector of length J containing the infection-age spline
#'   coefficients.
#' @param control An `rsv_control` object retained for interface consistency
#'   but not used directly in the current implementation.
#'
#' @return Numeric vector of length 365 containing the log visit-day masses.
#'
#' @details
#' For age day \eqn{m}, the function evaluates
#' \deqn{
#'   \log f_i^{\mathrm{visit}}(m)
#'   =
#'   \log c(m)
#'   +
#'   \log \lambda_i(m)
#'   +
#'   \log w(m)
#'   -
#'   H_i(m),
#' }
#' where
#' \deqn{
#'   H_i(m) = \beta^\top v_i(m).
#' }
#' The resulting event-time density values are used as visit-day masses in the
#' simulation distribution. Unlike the optimization objective for an observed
#' visit, this calculation retains the known \eqn{\log \lambda_i(m)} term.
logmass_i_visit <- function(lambda_i,
                            V_i,
                            w_vec,
                            c_vec,
                            beta,
                            control) {
  
  if (!is.numeric(lambda_i) || length(lambda_i) != 365L) {
    stop("lambda_i must be numeric vector length 365.")
  }
  
  if (!is.matrix(V_i) || ncol(V_i) != 366L) {
    stop("V_i must be numeric matrix with 366 columns.")
  }
  
  if (!is.numeric(beta) || length(beta) != nrow(V_i)) {
    stop("Length of beta must equal nrow(V_i).")
  }
  
  if (!is.numeric(w_vec) || length(w_vec) != 365L) {
    stop("w_vec must be numeric vector length 365.")
  }
  
  if (!is.numeric(c_vec) || length(c_vec) != 365L) {
    stop("c_vec must be numeric vector length 365.")
  }
  
  # Require finite positive values before applying logarithms.
  if (any(lambda_i <= 0 | !is.finite(lambda_i))) {
    stop("lambda_i must be strictly positive and finite.")
  }
  
  if (any(w_vec <= 0 | !is.finite(w_vec))) {
    stop("w_vec must be strictly positive and finite.")
  }
  
  if (any(c_vec <= 0 | !is.finite(c_vec))) {
    stop("c_vec must be strictly positive and finite.")
  }
  
  # Compute H_i(m) = beta^T v_i(m) for days 1 through 365.
  H_vec <- drop(crossprod(beta, V_i[, 2:366, drop = FALSE]))
  
  if (any(!is.finite(H_vec))) {
    stop("Non-finite integrated hazard values detected.")
  }
  
  # Evaluate the full log visit-density expression on the daily grid.
  log_mass_vec <-
    log(lambda_i) +
    log(c_vec) +
    log(w_vec) -
    H_vec
  
  log_mass_vec
}


# Visit outcome sampling ---------------------------------------------------


#' Sample healthcare-visit outcomes
#'
#' Samples one healthcare-visit outcome for each subject from a corresponding
#' subject-specific visit distribution and records the sampled visit age in
#' the RSV data object.
#'
#' @param data An `rsv_data` object containing the subjects whose
#'   healthcare-visit outcomes will be simulated.
#' @param visit_dist_list List of subject-specific visit distributions, with
#'   one numeric vector per subject. The first D entries represent visit days
#'   1 through D, and the final entry represents the no-visit outcome.
#' @param seed Integer scalar or `NULL`. If supplied, it is passed to
#'   `set.seed()` before sampling.
#' @param return_vectors Logical scalar indicating whether to return the
#'   sampled visit-age and visit-indicator vectors with the updated data.
#'
#' @return If `return_vectors = FALSE`, an `rsv_data` object with each
#'   subject's `visit_age` updated to the sampled visit day or `NA_integer_`
#'   for no visit.
#'
#'   If `return_vectors = TRUE`, a named list containing:
#' \describe{
#'   \item{\code{data}}{The updated `rsv_data` object.}
#'   \item{\code{visit_age}}{Integer vector containing the sampled visit ages,
#'     with `NA_integer_` for subjects with no visit.}
#'   \item{\code{I_visit}}{Integer vector containing 1 for subjects with a
#'     sampled visit and 0 for subjects with no visit.}
#' }
#'
#' @details
#' The number of modeled visit days D is inferred from the first distribution.
#' For each subject, a categorical outcome is sampled from the corresponding
#' distribution. Outcomes 1 through D are recorded as healthcare-visit ages,
#' while the final category is recorded as no visit.
#'
#' All subject-specific distributions must have the same length and contain
#' finite probability masses with positive total mass.
sample_visits_from_distribution <- function(data,
                                            visit_dist_list,
                                            seed = NULL,
                                            return_vectors = FALSE) {
  
  stopifnot(inherits(data, "rsv_data"))
  stopifnot(is.list(data$subjects), length(data$subjects) >= 1L)
  stopifnot(is.list(visit_dist_list))
  stopifnot(length(data$subjects) == length(visit_dist_list))
  stopifnot(is.logical(return_vectors), length(return_vectors) == 1L)
  
  n <- length(data$subjects)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Infer the number of visit days from the first distribution.
  first_dist <- visit_dist_list[[1]]
  stopifnot(is.numeric(first_dist), length(first_dist) >= 2L)
  
  # Reserve the final category for the no-visit outcome.
  D <- length(first_dist) - 1L  
  
  visit_age <- rep(NA_integer_, n)
  I_visit   <- integer(n)
  
  for (i in seq_len(n)) {
    p_i <- visit_dist_list[[i]]
    
    # Validate each supplied distribution before sampling.
    if (!is.numeric(p_i)) {
      stop(sprintf("visit_dist_list[[%d]] is not numeric.", i))
    }
    
    if (length(p_i) != D + 1L) {
      stop(sprintf(
        "visit_dist_list[[%d]] has length %d (expected %d).",
        i, length(p_i), D + 1L
      ))
    }
    
    if (any(!is.finite(p_i))) {
      stop(sprintf("Non-finite probabilities in visit_dist_list[[%d]].", i))
    }
    
    if (any(p_i < -1e-12)) {
      stop(sprintf("Negative probabilities in visit_dist_list[[%d]].", i))
    }
    
    s <- sum(p_i)
    if (!is.finite(s) || s <= 0) {
      stop(sprintf(
        "Invalid total mass in visit_dist_list[[%d]] (sum=%.16f).",
        i, s
      ))
    }
    
    # Sample one categorical outcome.
    k <- sample.int(D + 1L, size = 1L, prob = p_i)
    
    # Map visit-day outcomes to ages and the final category to no visit.
    if (k <= D) {
      visit_age[i] <- as.integer(k)
      I_visit[i]   <- 1L
    } else {
      visit_age[i] <- NA_integer_
      I_visit[i]   <- 0L
    }
    
    data$subjects[[i]]$visit_age <- visit_age[i]
  }
  
  # Preserve the rsv_data class explicitly.
  class(data) <- unique(c("rsv_data", class(data)))
  
  if (return_vectors) {
    return(list(
      data       = data,
      visit_age  = visit_age,
      I_visit    = I_visit
    ))
  }
  
  data
}
