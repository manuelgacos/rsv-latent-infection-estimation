# NOTE: Make this function more lightweight in memory for large simulations
# NOTE: Need to update the model (constant, periodic) functions to make the 
#   model check work
#' Generate subjects for RSV simulations (Bucket E)
#'
#' Creates a synthetic cohort of \code{n} subjects for the RSV simulation pipeline.
#' Each subject is assigned a calendar-time birthday index \code{birth_index} and
#' an initial \code{visit_age = NA_integer_} (to be filled later by the visit
#' simulation step).
#'
#' Birthdays are sampled uniformly on \eqn{\{1,\dots,\code{birth_max}\}}. To keep
#' the simulation components consistent, this function enforces compatibility
#' between \code{birth_max} and \code{model$age_len} when \code{age_len} is present
#' in the model.
#'
#' @param n Integer \eqn{\ge 1}. Number of subjects to generate.
#' @param model An \code{"rsv_model"} object. If \code{model$age_len} exists, it is
#'   used only for a consistency check against \code{birth_max}.
#' @param birth_max Integer > 0. Maximum birth day index (calendar days) used when
#'   sampling birthdays. This argument is required unless \code{model$age_len} is
#'   present and your implementation chooses to default \code{birth_max} from it.
#'
#'   \strong{Compatibility rule:}
#'   \itemize{
#'     \item If \code{model} does \emph{not} contain \code{age_len}, then
#'       \code{birth_max} must be provided (non-\code{NULL}).
#'     \item If \code{model} \emph{does} contain \code{age_len}, then
#'       \code{birth_max} must be provided and must satisfy
#'       \code{birth_max == model$age_len}; otherwise an error is thrown.
#'   }
#'
#' @param seed Optional integer seed for reproducibility. If supplied, the RNG
#'   state is restored on exit.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{data}}{An \code{"rsv_data"} object created by \code{make_rsv_data()}
#'     containing the generated subjects.}
#'   \item{\code{subjects}}{A list of subject records (as created by \code{make_subject()}).}
#'   \item{\code{birth_days}}{Integer vector of length \code{n} with sampled birthdays
#'     (the \code{birth_index} values).}
#' }
#'
#' @details
#' The generated subject records include:
#' \itemize{
#'   \item \code{id = 1:n}
#'   \item \code{birth_index} sampled uniformly from \code{1:birth_max}
#'   \item \code{visit_age = NA_integer_} (placeholder until Bucket F)
#' }
#'
#' @seealso \code{\link{make_sim_model_lambda_const}},
#'   \code{\link{make_sim_model_lambda_periodic}},
#'   \code{\link{simulate_visits_dataset}}
#'
#' @examples
#' \dontrun{
#' sim <- make_sim_model_lambda_const()
#' subj <- make_sim_subjects(n = 1000, model = sim$model, birth_max = 365, seed = 1)
#' head(subj$birth_days)
#' }
#'
#' @export
make_sim_subjects <- function(n, model, birth_max = NULL, seed = NULL) {
  # Need to update the model functions to make this work
  #stopifnot(inherits(model, "rsv_model"))
  stopifnot(is.numeric(n), length(n) == 1L, n >= 1)
  
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed)
  }
  
  # ---- Decide birth_max using "old vs new model" rules ----------------------
  has_age_len <- !is.null(model$age_len)
  
  if (!has_age_len) {
    # Old model: MUST provide birth_max
    if (is.null(birth_max)) {
      stop("make_sim_subjects: model has no 'age_len', so you must provide 'birth_max'.")
    }
  } else {
    # New model: MUST provide birth_max and it MUST match model$age_len
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
  
  # ---- Sample birth days (uniformly) and build subjects ---------------------------------
  birth_days <- sample.int(birth_max, size = as.integer(n), replace = TRUE)
  
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


#' Construct a single RSV subject object
#'
#' This function creates a standardized R list representing one subject in the
#' dataset. The resulting object contains all observable information required
#' for computing that subject's contribution to the log-likelihood.
#'
#' @param id Optional subject identifier (numeric, character, etc.).
#'
#' @param birth_index Integer. The subject's birth day expressed as a calendar
#'   index B_i. This index is used to shift the global lambda(t) curve onto the
#'   subject-specific age scale, i.e., lambda(a + B_i).
#'
#' @param visit_age Integer in 1:365, or NA. The age L_i (in days) at which the
#'   subject had their *first* bronchiolitis visit. If the subject had no visit
#'   in the first year of life, this should be NA. This value determines whether
#'   the subject belongs to I1 (visit before age 1) or I2 (no visit before age 1).
#'
#' @param weight Numeric scalar. Optional multiplicative weight for the
#'   log-likelihood (default = 1). This can be used to incorporate sampling
#'   weights or bootstrap weights if needed.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{id} – subject identifier.
#'     \item \code{birth_index} – integer birth index B_i.
#'     \item \code{visit_age} – integer L_i in 1:365, or NA if no visit by age 1.
#'     \item \code{weight} – weight applied to the subject’s likelihood term.
#'   }
#'
#' The returned object is meant to be stored inside a larger list of subjects
#' (via \code{make_rsv_data()}) and should be treated as read-only once created.
#'
make_subject <- function(id = NULL, birth_index, visit_age, weight = 1) {
  
  # ---- Input validation -----------------------------------------------------
  
  # birth_index must be a single finite number.
  # This index aligns the global lambda(t) with the subject's age scale.
  if (length(birth_index) != 1L || !is.finite(birth_index)) {
    stop("birth_index must be a single finite number.")
  }
  
  # visit_age must be either:
  #   * NA   (meaning no visit in the first year), OR
  #   * an integer in 1:365 (meaning the first visit occurred at that age).
  if (!is.na(visit_age)) {
    if (visit_age < 1 || visit_age > 365) {
      stop("visit_age must be in 1:365 or NA.")
    }
  }
  
  # ---- Return subject structure --------------------------------------------
  
  list(
    # Optional identifier, useful for debugging or reporting.
    id = id,
    
    # Birth day index B_i (stored as integer for consistency).
    birth_index = as.integer(birth_index),
    
    # First visit age L_i (NA_integer_ if no visit before age 1).
    visit_age = if (is.na(visit_age)) NA_integer_ else as.integer(visit_age),
    
    # Multiplicative weight for likelihood contribution.
    weight = weight
  )
}


#' Construct the RSV data object
#'
#' This function validates and packages a list of per-subject objects into a
#' standardized data container for use in the RSV likelihood. Each element of
#' the input \code{subjects} list must have been created by \code{make_subject()}
#' or must contain the same fields that \code{make_subject()} guarantees.
#'
#' @param subjects A list of subject objects, each created by
#'   \code{make_subject()}. Each subject must contain at least the fields:
#'   \itemize{
#'     \item \code{birth_index} – integer B_i.
#'     \item \code{visit_age}   – age at first visit L_i (1–365) or NA.
#'     \item \code{weight}      – likelihood weight.
#'   }
#'   Additional fields (e.g., \code{id}) are allowed.
#'
#' @return An object of class \code{"rsv_data"}:
#'   \itemize{
#'     \item \code{subjects} – the validated list of subject objects.
#'   }
#'
#' The resulting object is intended to be read-only in downstream computations.
#'
make_rsv_data <- function(subjects) {
  
  # ---- Basic structure checks ------------------------------------------------
  
  # 1. subjects must be a list
  if (!is.list(subjects)) {
    stop("`subjects` must be a list of subject objects.")
  }
  
  # Allow empty list (not typical, but consistent), but warn
  if (length(subjects) == 0L) {
    warning("`subjects` list is empty. Likelihood will be zero.")
  }
  
  # ---- Validate each subject -------------------------------------------------
  
  required_fields <- c("birth_index", "visit_age", "weight")
  
  for (i in seq_along(subjects)) {
    
    subj <- subjects[[i]]
    
    # Each subject must itself be a list
    if (!is.list(subj)) {
      stop(sprintf("subjects[[%d]] is not a list.", i))
    }
    
    # Check for required fields
    missing_fields <- setdiff(required_fields, names(subj))
    if (length(missing_fields) > 0L) {
      stop(sprintf(
        "subjects[[%d]] is missing required fields: %s",
        i,
        paste(missing_fields, collapse = ", ")
      ))
    }
    
    # No need to deeply validate values here (make_subject already handles that),
    # but we can enforce types lightly:
    if (!is.null(subj$birth_index) && !is.numeric(subj$birth_index)) {
      stop(sprintf("subjects[[%d]]$birth_index must be numeric/integer.", i))
    }
    if (!is.null(subj$weight) && !is.numeric(subj$weight)) {
      stop(sprintf("subjects[[%d]]$weight must be numeric.", i))
    }
  }
  
  # ---- Create the rsv_data object -------------------------------------------
  
  out <- list(subjects = subjects)
  class(out) <- "rsv_data"
  
  return(out)
}


#' Compute visit pmfs for all subjects (precomputed inputs)
#'
#' Returns a named list where each element is the full pmf
#' (length 366 named vector) for a subject.
#'
#' Assumes:
#'   - Global precompute (w_vec, c_vec) is already done
#'   - Subject precompute (lambda_i, V_i, id) is already done
#'
#' @param subj_pre_list list of subject precomputations
#'        Each element must contain:
#'          - id
#'          - lambda_i (length 365)
#'          - V_i (J x 366 matrix)
#' @param w_vec numeric vector length 365
#' @param c_vec numeric vector length 365
#' @param beta numeric vector length J
#' @param eta numeric vector length K
#' @param S_day numeric matrix K x 365
#' @param control rsv_control object
#'
#' @return named list of pmf vectors (one per subject)
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
  
  # ---- Assign names using subject IDs --------------------------------------
  
  # Requires rsv_precompute_subject() to be updated with IDs  
  # names(pmf_list) <- vapply(
  #   subj_pre_list,
  #   function(sp) as.character(sp$id),
  #   character(1)
  # )
  
  pmf_list
}


#' Full visit pmf for a single subject (likelihood-consistent)
#'
#' Returns a named probability vector of length 366:
#'   day001, ..., day365, no_visit
#'
#' The visit masses are computed using mass_i_visit(),
#' and the no-visit mass is computed using mass_i_no_visit(),
#' ensuring full consistency with the likelihood.
#'
#' @param lambda_i numeric vector length 365
#' @param V_i numeric matrix J x 366
#' @param w_vec numeric vector length 365
#' @param c_vec numeric vector length 365
#' @param beta numeric vector length J
#' @param eta numeric vector length K
#' @param S_day numeric matrix K x 365
#' @param control rsv_control object
#'
#' @return named numeric vector length 366 summing to 1
pmf_i <- function(lambda_i,
                  V_i,
                  w_vec,
                  c_vec,
                  beta,
                  eta,
                  S_day,
                  control, 
                  tol_prob = 1e-8) {
  
  # ---- Visit masses ---------------------------------------------------------
  
  visit_mass <- mass_i_visit(
    lambda_i = lambda_i,
    V_i      = V_i,
    w_vec    = w_vec,
    c_vec    = c_vec,
    beta     = beta,
    control  = control
  )
  
  # ---- No-visit mass --------------------------------------------------------
  
  tau_mass <- 1 - sum(visit_mass)
  
  # ---- Combine --------------------------------------------------------------
  
  pmf <- c(visit_mass, tau_mass)
  
  names(pmf) <- c(
    sprintf("day%03d", seq_len(365)),
    "no_visit"
  )
  
  # ---- Optional identity check ---------------------------------------------
  
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


#' Visit masses on natural scale (vectorized)
#'
#' Computes f_i(m) = lambda_i(m) * c(m) * w(m) * exp(-H_i(m))
#' for m = 1,...,365.
#'
#' This is a thin wrapper around logmass_i_visit() that
#' exponentiates the log-masses.
#'
#' @inheritParams logmass_i_visit
#'
#' @return numeric vector length 365 of visit masses
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


#   This functions disregard the functions from past implementations
#   including V2. 
# The main change in V4 was to compute the no visit probability
#   using the complement
#   NOTE: These functions also didn't work.
#   NOTE: The function worked with the w curve derived from Chris data.
#   NOTE: While testing, I discovered that using the functions in 
#     the likelihood don't yield a pmf since the sum is off 1. My
#     explanations lies within the implementation of the likelihood
#     in such a way that it is a dot product of beta and eta

#' Log visit-mass for all days (vectorized)
#'
#' Computes log f_i(m) for m = 1,...,365, where
#'
#'   f_i(m) = lambda_i(m) * c(m) * w(m) * exp(-H_i(m)),
#'
#' and H_i(m) = beta^T V_i[, m+1].
#'
#' This function is intended for pmf construction and restores
#' the lambda term omitted in loglik_i_visit() for optimization.
#'
#' @param lambda_i numeric vector length 365
#' @param V_i numeric matrix J x 366
#' @param w_vec numeric vector length 365
#' @param c_vec numeric vector length 365
#' @param beta numeric vector length J
#' @param control rsv_control object
#'
#' @return numeric vector length 365 of log visit-masses
logmass_i_visit <- function(lambda_i,
                            V_i,
                            w_vec,
                            c_vec,
                            beta,
                            control) {
  
  # ---- Basic shape checks (done once, not per-day) -------------------------
  
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
  
  # ---- Optional positivity check (cheap safety) ----------------------------
  
  if (any(lambda_i <= 0 | !is.finite(lambda_i))) {
    stop("lambda_i must be strictly positive and finite.")
  }
  
  if (any(w_vec <= 0 | !is.finite(w_vec))) {
    stop("w_vec must be strictly positive and finite.")
  }
  
  if (any(c_vec <= 0 | !is.finite(c_vec))) {
    stop("c_vec must be strictly positive and finite.")
  }
  
  # ---- Integrated hazard H_i(m) = beta^T V_i[, m+1] ------------------------
  
  H_vec <- drop(crossprod(beta, V_i[, 2:366, drop = FALSE]))
  
  if (any(!is.finite(H_vec))) {
    stop("Non-finite integrated hazard values detected.")
  }
  
  # ---- Log visit masses -----------------------------------------------------
  
  log_mass_vec <-
    log(lambda_i) +
    log(c_vec) +
    log(w_vec) -
    H_vec
  
  log_mass_vec
}


#' No-visit mass on natural scale (likelihood-consistent)
#'
#' Computes tau_i(beta, eta) on the natural scale using the exact
#' same computational pathway as loglik_i_no_visit().
#'
#' This is a thin wrapper around logmass_i_no_visit().
#'
#' @inheritParams logmass_i_no_visit
#'
#' @return numeric scalar tau_i
mass_i_no_visit <- function(lambda_i,
                            V_i,
                            w_vec,
                            beta,
                            eta,
                            S_day,
                            control) {
  
  log_tau <- logmass_i_no_visit(
    lambda_i = lambda_i,
    V_i      = V_i,
    w_vec    = w_vec,
    beta     = beta,
    eta      = eta,
    S_day    = S_day,
    control  = control
  )
  
  exp(log_tau)
}


#' Log no-visit mass (vectorized, likelihood-consistent)
#'
#' Computes log tau_i(beta, eta) using the exact same computational
#' pathway as loglik_i_no_visit(), but intended for pmf construction.
#'
#' @param lambda_i numeric vector length 365
#' @param V_i numeric matrix J x 366
#' @param w_vec numeric vector length 365
#' @param beta numeric vector length J
#' @param eta numeric vector length K
#' @param S_day numeric matrix K x 365
#' @param control rsv_control object
#'
#' @return numeric scalar log(tau_i)
logmass_i_no_visit <- function(lambda_i,
                               V_i,
                               w_vec,
                               beta,
                               eta,
                               S_day,
                               control) {
  
  # ---- Basic shape checks ---------------------------------------------------
  
  if (!is.numeric(lambda_i) || length(lambda_i) != 365L) {
    stop("lambda_i must be numeric vector length 365.")
  }
  
  if (!is.matrix(V_i) || !is.numeric(V_i) || ncol(V_i) != 366L) {
    stop("V_i must be numeric matrix with 366 columns.")
  }
  
  if (!is.numeric(beta) || length(beta) != nrow(V_i)) {
    stop("Length of beta must equal nrow(V_i).")
  }
  
  if (!is.numeric(w_vec) || length(w_vec) != 365L) {
    stop("w_vec must be numeric vector length 365.")
  }
  
  if (!is.numeric(eta)) {
    stop("eta must be numeric.")
  }
  
  if (!is.matrix(S_day) && !inherits(S_day, "Matrix")) {
    stop("S_day must be a numeric matrix (or Matrix).")
  }
  
  if (!is.numeric(S_day) || ncol(S_day) != 365L) {
    stop("S_day must be numeric matrix with 365 columns.")
  }
  
  if (length(eta) != nrow(S_day)) {
    stop("Length of eta must equal nrow(S_day).")
  }
  
  if (!inherits(control, "rsv_control") && !is.list(control)) {
    stop("control must be an 'rsv_control' object (or compatible list).")
  }
  
  check_bounds <- isTRUE(control$check_bounds)
  warn_on_clip <- isTRUE(control$warn_on_clip)
  eps_tau      <- if (!is.null(control$eps_tau)) control$eps_tau else 1e-12
  
  # ---- 1) Survival ----------------------------------------------------------
  
  Fbar <- Fbar_i(
    beta         = beta,
    V_i          = V_i,
    include_day0 = TRUE,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  if (!is.numeric(Fbar) || length(Fbar) != 366L) {
    stop("Fbar_i() must return numeric vector length 366.")
  }
  
  # ---- 2) Per-day infection probabilities ----------------------------------
  
  pi <- pi_from_w(
    lambda_shift_i = lambda_i,
    w              = w_vec,
    check_bounds   = check_bounds,
    warn_on_clip   = warn_on_clip
  )
  
  if (!is.numeric(pi) || length(pi) != 365L) {
    stop("pi_from_w() must return numeric vector length 365.")
  }
  
  # ---- 3) Visit kernel Q_i --------------------------------------------------
  
  Q <- Q_i(
    Fbar_i      = Fbar,
    pi_i        = pi,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  if (!is.numeric(Q) || length(Q) != 365L) {
    stop("Q_i() must return numeric vector length 365.")
  }
  
  # ---- 4) Aggregated kernel U_i --------------------------------------------
  
  U <- U_i(
    Q_i         = Q,
    S_day       = S_day,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  if (!is.numeric(U) || length(U) != length(eta)) {
    stop("U_i() must return numeric vector same length as eta.")
  }
  
  # ---- 5) tau_i and log -----------------------------------------------------
  
  tau <- tau_i(
    U_i         = U,
    eta         = eta,
    eps_tau     = eps_tau,
    check_bounds = check_bounds,
    warn_on_clip = warn_on_clip
  )
  
  if (!is.numeric(tau) || length(tau) != 1L ||
      !is.finite(tau) || tau <= 0) {
    stop("tau_i() must return a single positive finite scalar.")
  }
  
  log(tau)
}


sample_visits_from_distribution <- function(data,
                                            visit_dist_list,
                                            seed = NULL,
                                            return_vectors = FALSE) {
  # -----------------
  # Basic validation
  # -----------------
  stopifnot(inherits(data, "rsv_data"))
  stopifnot(is.list(data$subjects), length(data$subjects) >= 1L)
  stopifnot(is.list(visit_dist_list))
  stopifnot(length(data$subjects) == length(visit_dist_list))
  stopifnot(is.logical(return_vectors), length(return_vectors) == 1L)
  
  n <- length(data$subjects)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Infer number of days from first distribution
  first_dist <- visit_dist_list[[1]]
  stopifnot(is.numeric(first_dist), length(first_dist) >= 2L)
  
  D <- length(first_dist) - 1L  # last entry = no visit
  
  # Storage
  visit_age <- rep(NA_integer_, n)
  I_visit   <- integer(n)
  
  # -----------------
  # Main sampling loop
  # -----------------
  for (i in seq_len(n)) {
    p_i <- visit_dist_list[[i]]
    
    # ---- Defensive checks (no renormalization) ----
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
    
    # if (abs(s - 1) > 1e-10) {
    #   stop(sprintf(
    #     "visit_dist_list[[%d]] does not sum to 1 (sum=%.16f).",
    #     i, s
    #   ))
    # }
    
    # ---- Draw categorical sample ----
    k <- sample.int(D + 1L, size = 1L, prob = p_i)
    
    if (k <= D) {
      visit_age[i] <- as.integer(k)
      I_visit[i]   <- 1L
    } else {
      visit_age[i] <- NA_integer_
      I_visit[i]   <- 0L
    }
    
    # Update subject record
    data$subjects[[i]]$visit_age <- visit_age[i]
  }
  
  # Preserve class explicitly (defensive)
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