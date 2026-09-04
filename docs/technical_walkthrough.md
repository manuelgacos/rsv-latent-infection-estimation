# Technical Walkthrough

## Purpose and technical question

This walkthrough explains how the project estimates two hidden age-related patterns from partially observed healthcare data: the infection-age curve and the healthcare-visit weight.

The central question is:

> How can hidden first-infection and healthcare-visit patterns be estimated from birth timing, seasonal RSV circulation, and observed healthcare-visit outcomes?

Building on the broader infection-timing framework introduced by McKennan et al., this repository develops a day-level computational implementation for the partially observed setting. It connects the statistical model to the R code used for probability construction, penalized estimation, simulation, and numerical validation.

For a concise project overview, results, and instructions for running the example, see the main [README](../README.md).

## Observed and hidden processes

For subject $i$, let $B_i$ denote the birth index on the calendar-time grid. The known seasonal RSV circulation curve, $\lambda(t)$, describes how RSV activity changes over calendar time. A child's age $a$ corresponds to calendar time $B_i+a$, so children born on different dates may experience different circulation levels at the same age.

Let $R_i$ denote the child's age at first infection. This is the event the model seeks to understand, but it is not directly observed or stored in the estimation data.

Instead, the observed record contains $L_i$, the age of a healthcare visit during the first year of life, when one occurs. In the R implementation, `visit_age` is an integer from 1 through 365 for a child with a visit and `NA` for a child with no recorded visit during that period.

No visit does not imply no infection. A child may be infected without the infection producing a healthcare visit. The data therefore capture a selective consequence of the hidden process:

> Birth timing and seasonal circulation $\rightarrow$ hidden first infection $\rightarrow$ possible healthcare visit $\rightarrow$ observed visit age or no visit

The model must use these selectively observed outcomes to estimate the age-related component of first-infection risk and the healthcare-visit pattern.

### Key notation

| Symbol | Meaning |
| --- | --- |
| $B_i$                         | Birth index for subject $i$ on the calendar-time grid |
| $R_i$                         | Hidden age of first infection |
| $L_i$                         | Observed healthcare-visit age, when a visit occurs |
| $\lambda_i(a)=\lambda(B_i+a)$ | Seasonal RSV circulation experienced by subject $i$ at age $a$ |
| $w(a;\beta)$                  | Infection-age curve |
| $c(a;\eta)$                   | Nonnegative healthcare-visit weight |
| $Q_i(a;\beta)$                | Probability that subject $i$ is first infected during day $a$ |
| $\tau_i(\beta,\eta)$          | Probability that subject $i$ has no recorded visit during the first year |

## From birth timing to daily first-infection probabilities

The model is motivated in continuous time but evaluated on daily intervals. Day $a$ represents the interval $(a-1,a]$, for $a=1,\ldots,365$.

For a child with birth index $B_i$, the circulation level assigned to day $a$ is

$$
\lambda_i(a)=\lambda(B_i+a), \qquad a=1,\ldots,365.
$$

Birth timing is informative because children of the same age may encounter different circulation levels on different calendar dates.

The infection-age curve represents how susceptibility to first infection changes with age. The implementation begins with a continuous cubic B-spline basis and integrates each basis function over the corresponding one-day interval. The resulting vector $b(a)$ forms the stepwise basis for day $a$:

$$
w(a;\beta)=\beta^\top b(a),
$$

where $\beta$ is the vector of infection-age spline coefficients.

Within each day, seasonal circulation and the infection-age curve determine the conditional infection probability:

$$
\pi_i(a) = 1-\exp\left(-\lambda_i(a)w(a;\beta)\right).
$$

This is the probability of infection during day $a$, conditional on remaining uninfected through day $a-1$.

To calculate survival through age $d$, the implementation accumulates exposure over the preceding daily intervals. Let

$$
v_i(d) = \sum_{a=1}^{d}\lambda_i(a)\phi(a),
$$

where $\phi(a)$ is the day-integrated infection-age basis vector used in the cumulative calculation. In the implementation, $b(a)$ and $\phi(a)$ use the same daily integration rule but are stored separately for their different computational roles.

The cumulative infection hazard and survival probability are

$$
H_i(d;\beta)=\beta^\top v_i(d),
\qquad
\bar F_i(d;\beta)=\exp{-H_i(d;\beta)}.
$$

Thus, $\bar F_i(d;\beta)$ is the probability that subject $i$ remains uninfected through the end of day $d$.

The probability that the first infection occurs during day $a$ is

$$
Q_i(a;\beta) = \bar F_i(a-1;\beta)\pi_i(a;\beta).
$$

This combines survival through day $a-1$ with infection during day $a$. Across the first year, the values $Q_i(1;\beta),\ldots,Q_i(365;\beta)$ form the daily first-infection probability masses.

Daily quantities contain 365 values, while survival and cumulative-kernel objects also store their day-zero values, $\bar F_i(0;\beta)=1$ and $v_i(0)=0$, and therefore contain 366 entries.

## Connecting infection timing to healthcare visits

Infection time is hidden in the estimation data, which instead record either a healthcare visit age $L_i$ or no visit during the first year.

The model connects infection timing to the observed outcomes through

$$
c(a;\eta)=\eta^\top s(a),
$$

where $s(a)$ is the day-level healthcare-visit basis vector and $\eta$ is the vector of healthcare-visit spline coefficients.

The curve $c(a;\eta)$ is a nonnegative healthcare-visit weight. It allows infections at different ages to contribute differently to the observed visit process, but it is not a standalone probability bounded between zero and one.

For a subject with a recorded visit at age $L_i$, the continuous-time model gives the visit density

$$
f_i^{\mathrm{visit}}(L_i;\beta,\eta) 
= c(L_i;\eta)\lambda_i(L_i)w(L_i;\beta) \exp{-H_i(L_i;\beta)}.
$$

The corresponding full log-likelihood contribution is

$$
\ell_i^{\mathrm{visit}}(\beta,\eta)
= \log c(L_i;\eta) + \log \lambda_i(L_i) + \log w(L_i;\beta) - H_i(L_i;\beta).
$$

The implementation evaluates this density at the recorded visit day. It does not include a separate delay distribution between infection and the associated healthcare visit. Because $\lambda_i(L_i)$ is known and does not depend on either coefficient vector, its logarithm is omitted from the computational optimization objective.

For a subject with no recorded visit, the contribution is the probability that no healthcare visit occurs during the first year. Under the daily stepwise representation, the implementation aggregates the possible first-infection days:

$$
U_i(\beta) = \sum_{a=1}^{365}s(a)Q_i(a;\beta).
$$

The resulting first-year visit probability is

$$
g_i(\beta,\eta) 
= \eta^\top U_i(\beta) 
= \sum_{a=1}^{365}c(a;\eta)Q_i(a;\beta),
$$

and its complement is the no-visit probability:

$$
\tau_i(\beta,\eta) = 1-g_i(\beta,\eta).
$$

The corresponding log-likelihood contribution is

$$
\ell_i^{\mathrm{no\ visit}}(\beta,\eta) = \log\tau_i(\beta,\eta).
$$

With subject-level weight $\omega_i$, the total log-likelihood is

$$
\ell(\beta,\eta) = \sum_{i=1}^{n}\omega_i\ell_i(\beta,\eta).
$$

The reference example assigns every subject a weight of one.

Observed visits therefore contribute event-time densities, whereas no-visit records contribute first-year no-visit probabilities. In simulation, the visit density is evaluated on the daily grid to construct visit-day masses and the complementary no-visit outcome.

## Representing and regularizing the two curves

The day-integrated spline bases allow both curves to vary over age without estimating an unrelated value for every day:

$$
w(a;\beta)=\beta^\top b(a),
\qquad
c(a;\eta)=\eta^\top s(a).
$$

The reference example uses 14 basis functions for each curve.

Because the day-level basis values are nonnegative, constraining every coefficient to be nonnegative keeps both estimated curves nonnegative across the age grid:

$$
\beta \geq 0,
\qquad
\eta \geq 0.
$$

These constraints control the curves themselves. They do not alone guarantee that every combined subject-level visit probability is valid, so those probabilities are also checked and stabilized during likelihood evaluation.

To discourage unnecessarily sharp changes, the implementation penalizes second differences between adjacent coefficients. For a coefficient vector $\theta$, a second difference is

$$
\theta_{j+2}-2\theta_{j+1}+\theta_j.
$$

Large second differences indicate greater local variation in the coefficient sequence.

Let $M_\beta$ and $M_\eta$ denote the second-difference penalty matrices for the two coefficient vectors, and let $\alpha \geq 0$ denote the common smoothing parameter. The full penalized objective is

$$
F(\beta,\eta)
= -\ell(\beta,\eta) + \frac{\alpha}{2} 
\left( \beta^\top M_\beta\beta + \eta^\top M_\eta\eta \right),
$$

subject to $\beta \geq 0$ and $\eta \geq 0$.

The negative log-likelihood measures disagreement with the observed outcomes, while the quadratic terms penalize roughness. Larger values of $\alpha$ favor smoother curves, while smaller values allow greater local variation.

The reference example fixes $\alpha=50$. Selecting the smoothing level through cross-validation is outside the scope of this repository.

## Analytic gradients and alternating estimation

Because $\beta$ and $\eta$ are connected through the likelihood, the implementation uses analytic gradients to reduce the cost of repeated optimization across subjects and possible infection days.

For a subject with an observed visit, the gradients of the log-likelihood are

$$
\nabla_\beta \ell_i^{\mathrm{visit}} = \frac{b(L_i)}{w(L_i;\beta)} - v_i(L_i),
\qquad
\nabla_\eta \ell_i^{\mathrm{visit}} = \frac{s(L_i)}{c(L_i;\eta)}.
$$

The $\eta$ gradient depends on the healthcare-visit basis and curve value at $L_i$, while the $\beta$ gradient also accounts for cumulative infection exposure through that day.

For a subject with no observed visit, the gradient with respect to $\eta$ is

$$
\nabla_\eta \ell_i^{\mathrm{no\ visit}} = -\frac{U_i(\beta)}{\tau_i(\beta,\eta)}.
$$

The gradient with respect to $\beta$ is

$$
\nabla_\beta \ell_i^{\mathrm{no\ visit}} = -\frac{1}{\tau_i(\beta,\eta)}
\sum_{a=1}^{365} c(a;\eta) \left[ \bar F_i(a;\beta)\lambda_i(a)b(a) - Q_i(a;\beta)v_i(a-1) \right].
$$

This expression accounts for how changing $\beta$ affects both survival and the probability assigned to each possible first-infection day. The implementation evaluates its main terms through matrix-vector operations rather than an explicit loop over all 365 ages.

The subject-level gradients are combined using the subject weights. The implementation then negates the resulting log-likelihood gradient and adds the penalty gradients to obtain the gradient of the minimized objective.

Estimation alternates between the two parameter blocks:

1. Update $\beta$ while holding $\eta$ fixed.
2. Update $\eta$ while holding the updated $\beta$ fixed.
3. Evaluate the full penalized objective and repeat.

Both updates use L-BFGS-B, a bounded numerical optimizer, with lower bounds of zero. This alternating procedure is a form of block coordinate descent.

The outer loop supports objective-based, coefficient-based, or combined convergence criteria. The reference example stops when the relative objective change falls below $10^{-5}$.

## Computational design and numerical safeguards

Repeated likelihood and gradient evaluations are expensive because the model performs daily probability calculations for every subject. The implementation reduces this cost by precomputing quantities that do not change with $\beta$ or $\eta$.

For each subject, the code calculates once and reuses:

* The subject-specific circulation curve $\lambda_i(a)$.
* The cumulative kernel $v_i(d)$ for days 0 through 365.

These quantities depend on the subject's birth index, the fixed seasonal circulation curve, and the spline basis. During estimation, only the parameter-dependent curves, survival probabilities, infection-day masses, and likelihood contributions must be recomputed.

Numerical reliability is handled separately. Quantities that are theoretically constrained are checked and stabilized when necessary. Small negative values caused by floating-point error can be clipped to zero, probability-like quantities can be restricted to valid ranges, and very small positive values are bounded away from zero before logarithms are evaluated. The no-visit probability $\tau_i(\beta,\eta)$ is also assigned a positive lower bound to prevent undefined or unstable log-likelihood values.

Controls such as `check_bounds` and `warn_on_clip` can be configured independently. Bound checking determines whether substantial violations cause errors, while warning controls determine whether clipping is reported. The reference example sets both flags to `FALSE` and validates the resulting subject-level distributions separately.

The workflow verifies that the simulated distributions are finite, nonnegative, and normalized within numerical tolerance. During estimation, it records whether the inner optimizers converge and terminates the alternating procedure if the full objective becomes non-finite.

## Simulation and recovery validation

Because first-infection ages are hidden in the estimation data, the reference example uses simulation to evaluate whether the workflow can recover known age-related patterns from visit and no-visit outcomes.

The example loads a fixed seasonal RSV circulation curve and defines target infection-age and healthcare-visit curves. These targets are projected onto the model's spline bases so that data generation and estimation use the same day-level representation.

The simulation then:

1. Generates 20,000 subjects with different birth indices.
2. Aligns each subject's first year with the seasonal circulation curve.
3. Evaluates the visit density at the 365 daily grid points.
4. Uses these values as visit-day masses and assigns the remaining mass to the no-visit outcome.
5. Samples one observed outcome for each subject.
6. Estimates $\beta$ and $\eta$ from the birth indices, circulation curve, and observed visit records.

Before sampling, each resulting 366-entry distribution is checked for finite values, nonnegative masses, and normalization.

The known simulation coefficients are perturbed using a fixed random seed to create reproducible starting values. They are not supplied to the optimizer as targets to match. After estimation, the known curves are used only to evaluate recovery.

The reference run produced 5,989 observed visits, corresponding to a visit rate of 29.9%. The mean model-implied probability of first infection by age one was 50.3%. The alternating procedure satisfied its objective-based convergence criterion after three outer iterations, and the complete workflow took approximately seven to eight minutes.

| Validation measure      | Infection-age curve | Healthcare-visit curve |
| ----------------------- | ------------------: | ---------------------: |
| Root mean squared error |            0.000453 |                0.01184 |
| Maximum absolute error  |            0.000891 |                0.02366 |
| Correlation             |              0.9914 |                 0.9974 |

Root mean squared error summarizes the typical daily discrepancy between an estimated curve and its simulation truth. Maximum absolute error records the largest local difference, while correlation measures agreement in overall shape. Because the curves operate on different numerical scales, their error magnitudes should be interpreted within each curve rather than compared directly.

The maximum error in the sums of the subject-level distributions was $1.11\times10^{-16}$, confirming normalization to one within floating-point precision.

Together, these results show that the reference configuration recovers the broad shapes of both hidden curves.

## Mapping the method to the repository

The repository modules follow the same sequence as the statistical workflow.

* **`R/model_setup.R`** constructs the spline bases, penalty matrices, model objects, and numerical-control settings.

* **`R/probability_core.R`** aligns circulation with each subject's age and evaluates cumulative exposure, survival, daily infection probabilities, infection-day masses, and visit probabilities. `rsv_precompute_subject()` prepares fixed subject-level quantities for reuse during estimation.

* **`R/likelihood.R`** implements the observed-visit density and no-visit probability contributions.

* **`R/estimation.R`** contains the analytic gradients, penalized objectives, nonnegative blockwise optimizers, and `estimate_model_alternating()`.

* **`R/simulation.R`** generates subjects, constructs visit and no-visit distributions, and samples the observed outcomes.

* **`R/diagnostics.R`** provides helpers for summarizing outcomes, calculating probabilities, scaling reference curves, and creating recovery plots.

* **`example/run_example.R`** coordinates the complete reference workflow. It loads the inputs, simulates the data, validates the subject-level distributions, estimates the model, calculates recovery metrics, reports convergence and runtime, and saves the comparison figures.

## Scope and technical limitations

The reference example uses simulated data and one fixed configuration so that the hidden infection-age and healthcare-visit curves are known and recovery can be measured directly. The smoothing parameter is fixed, and the public workflow does not include cross-validation, uncertainty intervals, or evaluation across multiple sample sizes, visit rates, seasonal patterns, or curve shapes.

The reported results therefore demonstrate recovery for the included configuration rather than guaranteeing comparable behavior in other settings. Applying the framework to observational healthcare data would also require data preparation, definitions of eligible visits and follow-up windows, validation of the circulation input, assessment of model assumptions, and investigation of selection and measurement bias.

The repository is a focused implementation of the core simulation, estimation, and validation workflow rather than the complete research codebase or a general-purpose R package. Within this scope, it demonstrates the complete path from a partially observed statistical model to reproducible simulation, constrained estimation, and numerical recovery validation.