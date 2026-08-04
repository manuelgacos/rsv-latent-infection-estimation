# Reproducible simulation and estimation example for the RSV model.
# Run this script from the root of the R project.

total_start <- Sys.time()


# 01. Setup ---------------------------------------------------------------

# Check dependencies and load the project functions used in this example.

required_packages <- c("splines2")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following package(s) before running the example: ",
    paste(missing_packages, collapse = ", ")
  )
}

source_files <- c(
  file.path("R", "model_setup.R"),
  file.path("R", "probability_core.R"),
  file.path("R", "likelihood.R"),
  file.path("R", "estimation.R"),
  file.path("R", "simulation.R"),
  file.path("R", "diagnostics.R")
)

missing_source_files <- source_files[!file.exists(source_files)]

if (length(missing_source_files) > 0L) {
  stop(
    paste0(
      "The expected project files were not found.\n\n",
      "Run this script from the project root.\n\n",
      "Current working directory:\n  ",
      getwd(),
      "\n\nMissing files:\n",
      paste0("  - ", missing_source_files, collapse = "\n")
    )
  )
}

library(splines2)

invisible(lapply(source_files, source))


# 02. Configuration -------------------------------------------------------

# Define the reproducible simulation and estimation settings.

config <- list(
  # Random seeds
  seed_subjects = 1153318542L,
  seed_visits   = 1485562642L,
  seed_init     = 1354052710L,
  
  # Simulation size
  n_subjects = 20000L,
  
  # Time grids
  age_days     = 365L,
  birth_max    = 578L,
  calendar_len = 943L,
  
  # Spline basis
  J      = 14L,
  K      = 14L,
  degree = 3L,
  
  # Estimation
  alpha     = 50,
  noise_sd  = 0.05,
  max_outer = 10L,
  tol_obj   = 1e-5,
  tol_l2    = 1e-6,
  
  # Output
  save_plots = TRUE
)


# 03. Load fixed inputs ---------------------------------------------------

# Load the fixed circulation curve, reference curve, and scaling constants.

input_files <- c(
  lambda_global = file.path("data", "lambda_global.rds"),
  w_reference   = file.path("data", "w_reference.rds"),
  w_scale       = file.path("data", "w_scale.rds"),
  c_scale       = file.path("data", "c_scale.rds")
)

missing_inputs <- input_files[!file.exists(input_files)]

if (length(missing_inputs) > 0L) {
  stop(
    paste0(
      "The required input files were not found:\n",
      paste0("  - ", missing_inputs, collapse = "\n")
    )
  )
}

lambda_global <- readRDS(input_files[["lambda_global"]])
w_reference   <- readRDS(input_files[["w_reference"]])
w_scale       <- readRDS(input_files[["w_scale"]])
c_scale       <- readRDS(input_files[["c_scale"]])

if (
  !is.numeric(lambda_global) ||
  length(lambda_global) != config$calendar_len ||
  any(!is.finite(lambda_global)) ||
  any(lambda_global < 0)
) {
  stop(
    "`lambda_global` must be a nonnegative numeric vector of length ",
    config$calendar_len,
    "."
  )
}

if (
  !is.numeric(w_reference) ||
  length(w_reference) != config$age_days ||
  any(!is.finite(w_reference)) ||
  any(w_reference < 0)
) {
  stop(
    "`w_reference` must be a nonnegative numeric vector of length ",
    config$age_days,
    "."
  )
}

if (
  !is.numeric(w_scale) ||
  length(w_scale) != 1L ||
  !is.finite(w_scale) ||
  w_scale <= 0
) {
  stop("`w_scale` must be a single positive finite numeric value.")
}

if (
  !is.numeric(c_scale) ||
  length(c_scale) != 1L ||
  !is.finite(c_scale) ||
  c_scale <= 0
) {
  stop("`c_scale` must be a single positive finite numeric value.")
}


# 04. Construct the model -------------------------------------------------

# Configure numerical checks and warning behavior for the example.
control <- make_rsv_control(
  check_bounds = FALSE,
  warn_on_clip = FALSE
)

# Build the spline bases and combine them with the fixed circulation curve.
model_setup <- make_sim_model_lambda_custom(
  lambda_global = lambda_global,
  J              = config$J,
  K              = config$K,
  degree         = config$degree,
  days           = config$age_days,
  age_len        = config$age_days,
  birth_max      = config$birth_max,
  calendar_len   = config$calendar_len
)

model <- model_setup$model

if (!inherits(model, "rsv_model")) {
  stop("Model construction did not return a valid `rsv_model` object.")
}


# 05. Construct true curves -----------------------------------------------

# Scale the reference infection-age curve to the target infection level.
w_target <- w_reference * w_scale

# Define a decreasing visit-age curve and scale it to the target visit level.
c_target <- seq(
  from       = 1,
  to         = 0.5,
  length.out = config$age_days
) * c_scale

# Project both target curves onto the model's spline bases.
truth <- fit_true_curves(
  model     = model,
  w_target  = w_target,
  c_target  = c_target,
  days      = config$age_days,
  normalize = FALSE,
  opt_maxit = 2000L
)

beta_true <- truth$beta_true
eta_true  <- truth$eta_true

w_true <- truth$w_true
c_true <- truth$c_true

if (
  any(!is.finite(c(beta_true, eta_true, w_true, c_true))) ||
  any(beta_true < 0) ||
  any(eta_true < 0)
) {
  stop("The true curves could not be represented by valid spline coefficients.")
}


# 06. Generate subjects ---------------------------------------------------

# Simulate birth dates before healthcare visits are observed.
subject_setup <- make_sim_subjects(
  n         = config$n_subjects,
  model     = model,
  birth_max = config$birth_max,
  seed      = config$seed_subjects
)

data_unobserved <- subject_setup$data
subjects <- subject_setup$subjects

if (
  !inherits(data_unobserved, "rsv_data") ||
  length(subjects) != config$n_subjects
) {
  stop("Subject generation did not return the expected simulated cohort.")
}


# 07. Precompute subject quantities --------------------------------------

# Precompute subject-specific quantities reused in simulation and estimation.
precompute_start <- Sys.time()

subj_pre_list <- lapply(
  subjects,
  rsv_precompute_subject,
  model   = model,
  control = control
)

precompute_end <- Sys.time()

precompute_runtime <- difftime(
  precompute_end,
  precompute_start,
  units = "mins"
)

if (length(subj_pre_list) != config$n_subjects) {
  stop("Subject-level precomputation did not return one result per subject.")
}


# 08. Construct visit distributions --------------------------------------

# Compute the true infection-age and visit-age curves on the daily grid.
global_truth <- rsv_precompute_global(
  beta    = beta_true,
  eta     = eta_true,
  model   = model,
  control = control
)

# Construct each subject's probability distribution for visit age or no visit.
pmf_tolerance <- 1e-8

visit_distributions <- pmf_dataset(
  subj_pre_list = subj_pre_list,
  w_vec         = global_truth$w_vec,
  c_vec         = global_truth$c_vec,
  beta          = beta_true,
  eta           = eta_true,
  S_day         = model$S_day,
  control       = control,
  tol_prob      = pmf_tolerance
)

pmf_sums <- vapply(
  visit_distributions,
  sum,
  numeric(1)
)

max_pmf_error <- max(abs(pmf_sums - 1))

if (
  length(visit_distributions) != config$n_subjects ||
  any(!is.finite(pmf_sums)) ||
  max_pmf_error > pmf_tolerance
) {
  stop("The subject-level visit distributions are not valid probability distributions.")
}


# 09. Simulate healthcare visits -----------------------------------------

# Sample each subject's observed visit outcome from the visit distributions.
data_observed <- sample_visits_from_distribution(
  data            = data_unobserved,
  visit_dist_list = visit_distributions,
  seed            = config$seed_visits,
  return_vectors  = FALSE
)

if (
  !inherits(data_observed, "rsv_data") ||
  length(data_observed$subjects) != config$n_subjects
) {
  stop("Visit simulation did not return the expected observed dataset.")
}


# 10. Validate the simulated data ----------------------------------------

# Summarize the observed visits and model-implied infection probabilities.
visit_summary <- extract_visit_ages(
  data          = data_observed,
  model         = model,
  visit_age_max = config$age_days,
  days          = config$age_days,
  compute_g     = FALSE
)

infection_probabilities <- vapply(
  subj_pre_list,
  prob_infected_by_age,
  numeric(1),
  beta    = beta_true,
  days    = config$age_days,
  control = control
)

simulation_summary <- list(
  n_subjects = config$n_subjects,
  n_visits = sum(visit_summary$I_visit),
  visit_rate = mean(visit_summary$I_visit),
  mean_infection_probability = mean(infection_probabilities),
  max_pmf_error = max_pmf_error
)

cat(
  "\nSimulation summary\n",
  "------------------\n",
  sprintf("%-39s %s\n",
          "Subjects:",
          format(simulation_summary$n_subjects, big.mark = ",")),
  sprintf("%-39s %s\n",
          "Observed visits:",
          format(simulation_summary$n_visits, big.mark = ",")),
  sprintf("%-39s %.1f%%\n",
          "Observed visit rate:",
          100 * simulation_summary$visit_rate),
  sprintf("%-39s %.1f%%\n",
          "Mean infection probability by age 1:",
          100 * simulation_summary$mean_infection_probability),
  sprintf("%-39s %.2e\n",
          "Maximum PMF error:",
          simulation_summary$max_pmf_error),
  sep = ""
)


# 11. Fit the model -------------------------------------------------------

# Create reproducible starting values by perturbing the true coefficients.
initial_values <- make_perturbed_initial_values(
  beta     = beta_true,
  eta      = eta_true,
  noise_sd = config$noise_sd,
  seed     = config$seed_init
)

beta_init <- initial_values$beta_init
eta_init  <- initial_values$eta_init

cat(
  "\nModel estimation\n",
  "----------------\n",
  sep = ""
)

# Estimate the infection-age and visit-age curves by alternating optimization.
fit_start <- Sys.time()

fit <- estimate_model_alternating(
  beta_init     = beta_init,
  eta_init      = eta_init,
  alpha         = config$alpha,
  data          = data_observed,
  model         = model,
  rsv_control   = control,
  subj_pre_list = subj_pre_list,
  max_outer     = config$max_outer,
  convergence   = "objective",
  tol_obj       = config$tol_obj,
  tol_l2        = config$tol_l2,
  verbose       = TRUE
)

fit_end <- Sys.time()

fit_runtime <- difftime(
  fit_end,
  fit_start,
  units = "mins"
)

if (
  !is.list(fit) ||
  any(!is.finite(c(fit$beta_hat, fit$eta_hat, fit$final_objective)))
) {
  stop("Model estimation did not return usable results.")
}

if (!isTRUE(fit$converged)) {
  warning(
    "The alternating optimization stopped without satisfying the ",
    "convergence criterion after ",
    fit$n_outer,
    " outer iteration(s)."
  )
}


# 12. Extract and evaluate estimates -------------------------------------

# Evaluate the fitted infection-age and visit-age curves on the daily grid.
global_estimate <- rsv_precompute_global(
  beta    = fit$beta_hat,
  eta     = fit$eta_hat,
  model   = model,
  control = control
)

w_hat <- global_estimate$w_vec
c_hat <- global_estimate$c_vec

if (
  any(!is.finite(c(w_hat, c_hat))) ||
  length(w_hat) != config$age_days ||
  length(c_hat) != config$age_days
) {
  stop("The fitted coefficients did not produce valid estimated curves.")
}

# Compare the estimated curves with the simulation truth.
curve_metrics <- list(
  w_rmse        = sqrt(mean((w_hat - w_true)^2)),
  c_rmse        = sqrt(mean((c_hat - c_true)^2)),
  w_max_error   = max(abs(w_hat - w_true)),
  c_max_error   = max(abs(c_hat - c_true)),
  w_correlation = cor(w_hat, w_true),
  c_correlation = cor(c_hat, c_true)
)

cat(
  "\nCurve recovery\n",
  "--------------\n",
  sprintf("%-30s %.4e\n",
          "Infection curve RMSE:",
          curve_metrics$w_rmse),
  sprintf("%-30s %.4e\n",
          "Infection maximum error:",
          curve_metrics$w_max_error),
  sprintf("%-30s %.4f\n",
          "Infection correlation:",
          curve_metrics$w_correlation),
  sprintf("%-30s %.4e\n",
          "Visit curve RMSE:",
          curve_metrics$c_rmse),
  sprintf("%-30s %.4e\n",
          "Visit maximum error:",
          curve_metrics$c_max_error),
  sprintf("%-30s %.4f\n",
          "Visit correlation:",
          curve_metrics$c_correlation),
  sep = ""
)

curve_results <- data.frame(
  age_days = seq_len(config$age_days),
  age_months = seq_len(config$age_days) / (365 / 12),
  w_true = w_true,
  w_estimated = w_hat,
  c_true = c_true,
  c_estimated = c_hat
)


# 13. Plot results and summarize runtime ---------------------------------

# Plot the true and estimated curves.
plot_curve_recovery(
  curve_results = curve_results,
  save_plots    = config$save_plots,
  output_dir    = "outputs"
)

# Record the end of the complete workflow.
total_end <- Sys.time()

total_runtime <- difftime(
  total_end,
  total_start,
  units = "mins"
)

runtime_summary <- list(
  precomputation_minutes = as.numeric(precompute_runtime),
  estimation_minutes     = as.numeric(fit_runtime),
  total_minutes          = as.numeric(total_runtime)
)

cat(
  "\nRuntime summary\n",
  "---------------\n",
  sprintf(
    "%-25s %.2f minutes\n",
    "Subject precomputation:",
    runtime_summary$precomputation_minutes
  ),
  sprintf(
    "%-25s %.2f minutes\n",
    "Model estimation:",
    runtime_summary$estimation_minutes
  ),
  sprintf(
    "%-25s %.2f minutes\n",
    "Total runtime:",
    runtime_summary$total_minutes
  ),
  sep = ""
)
