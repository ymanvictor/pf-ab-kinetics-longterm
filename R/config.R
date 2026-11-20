# Settings shared by all scripts: antigens and antibodies, which combinations
# to fit, the simulated study design, priors, sampler settings and the
# definitions used for seropositivity and seroreversion.
#
# Time unit in the models: months of 30 days. Years in the outputs are
# 12 model months.

# ------------------------------------------------------------------------------
# Antigens, antibodies and models
# ------------------------------------------------------------------------------

ANTIGENS   <- c("MSP1", "MSP2", "MSP3", "AMA1", "GAMA")
ANTIBODIES <- c("IgM", "IgG", "IgG1", "IgG3")

# IgM is fitted with the two-component model, IgG and its subclasses with the
# three-component model.
model_for_antibody <- function(antibody) {
  ifelse(antibody == "IgM", "2-component", "3-component")
}

STAN_FILES <- c(
  "2-component" = "stan/Ab_kin_2comp.stan",
  "3-component" = "stan/Ab_kin_3comp.stan"
)

# Combinations fitted by 02_fit_models.R and summarised by 03 and 04.
# Data are simulated for all 20 combinations; by default two are fitted, one
# per model. Set FIT_ALL_COMBINATIONS <- TRUE to fit all 20.
FIT_ALL_COMBINATIONS <- FALSE

COMBINATIONS_TO_FIT <- data.frame(
  antigen  = c("MSP1", "MSP1"),
  antibody = c("IgM",  "IgG")
)

if (FIT_ALL_COMBINATIONS) {
  COMBINATIONS_TO_FIT <- expand.grid(
    antigen = ANTIGENS,
    antibody = ANTIBODIES,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Simulated study design (01_simulate_data.R)
# ------------------------------------------------------------------------------
# Mimics the structure of the study, not its actual sampling times:
#   - cases sampled densely from presentation over the first year, plus one
#     sample 5-10 years later;
#   - cases seen only once to three times, 8-45 years after infection;
#   - unexposed controls with a single sample each.

SIM_SEED <- 2026

N_DENSE_CASES <- 14
N_LATE_CASES  <- 32
N_CONTROLS    <- 30

DENSE_VISIT_DAYS      <- c(0, 7, 30, 90, 180, 365)  # nominal visits
DENSE_VISIT_JITTER    <- 0.2                        # +/- 20% of the visit day
P_MISSED_VISIT        <- 0.08                       # each visit after day 0
DENSE_LATE_YEARS      <- c(5, 10)                   # one late sample
LATE_FIRST_YEARS      <- c(8, 45)                   # late-only cases
LATE_N_SAMPLES_PROBS  <- c(0.5, 0.35, 0.15)         # 1, 2 or 3 samples
LATE_GAP_YEARS        <- c(1, 10)                   # between repeat samples
MAX_FOLLOW_UP_YEARS   <- 45

DAYS_PER_YEAR <- 365.25

# Upper limit of quantification on the nMFI scale. In the study the limit was
# a raw MFI of 28,000, converted to nMFI per plate; here it is one value.
UPPER_LIMIT_NMFI <- 28000

# Delays (months) are drawn from the prior for densely sampled cases and are
# zero otherwise, as in the model.
DELTA_SD    <- 1 / 3
DELTA_RANGE <- c(-1, 1)

# Population parameters used to simulate the data: posterior medians from the
# paper's Supplementary Data 1.
SIM_PARAMETER_FILE <- "data/Supplementary_Data_1_posterior_parameter_summary.csv"

# ------------------------------------------------------------------------------
# Data handling in the fits (02_fit_models.R)
# ------------------------------------------------------------------------------

DAYS_PER_MONTH <- 30

# Delays are estimated only for cases with at least DELTA_MIN_EARLY_OBS
# samples within DELTA_ESTIMATE_CUTOFF_MONTHS of presentation.
DELTA_ESTIMATE_CUTOFF_MONTHS <- 3
DELTA_MIN_EARLY_OBS          <- 2L

# Early-peak window for the data-based log_gB_pop prior (peak rule).
GB_PRIOR_PEAK_WINDOW_MONTHS <- 3
LOG_GB_PRIOR_SD             <- 1

# ------------------------------------------------------------------------------
# Priors
# ------------------------------------------------------------------------------
# Rates are per month. Kinetic priors are normal on the log rate (centred on a
# half-life) or on the logit weight.

log_rate_from_half_life_days <- function(days) {
  log(DAYS_PER_MONTH * log(2) / days)
}

log_rate_from_daily_rate <- function(rate_per_day) {
  log(DAYS_PER_MONTH * rate_per_day)
}

KINETIC_PRIORS <- list(
  IgM = list(
    log_r_prior_loc       = log_rate_from_half_life_days(5),
    log_r_prior_sd        = 0.15,
    log_c_long_prior_loc  = log_rate_from_half_life_days(45),
    log_c_long_prior_sd   = 0.65,
    log_c_short_prior_loc = log_rate_from_half_life_days(10),
    log_c_short_prior_sd  = 0.50,
    rho_prior_loc         = qlogis(0.80),
    rho_prior_sd          = 0.80
  ),
  IgG = list(
    log_r_prior_loc              = log_rate_from_half_life_days(21),
    log_r_prior_sd               = 0.15,
    log_c_long_prior_loc         = log_rate_from_half_life_days(1825),
    log_c_long_prior_sd          = 0.60,
    log_c_intermediate_prior_loc = log_rate_from_daily_rate(0.0025),
    log_c_intermediate_prior_sd  = 0.60,
    log_c_short_prior_loc        = log_rate_from_half_life_days(10),
    log_c_short_prior_sd         = 0.50,
    rho_prior_loc                = qlogis(0.80),
    rho_prior_sd                 = 1.50,
    gamma_prior_loc              = qlogis(0.60),
    gamma_prior_sd               = 0.70
  )
)
KINETIC_PRIORS$IgG1 <- KINETIC_PRIORS$IgG
KINETIC_PRIORS$IgG3 <- modifyList(
  KINETIC_PRIORS$IgG,
  list(log_r_prior_loc = log_rate_from_half_life_days(7))
)

# Half-normal scales for the between-individual SDs and the residual SD.
RE_SD_PRIORS <- list(
  "2-component" = list(
    A0_sigma_prior_sd      = 1.00,
    gB_sigma_prior_sd      = 2.00,
    c_long_sigma_prior_sd  = 0.35,
    c_short_sigma_prior_sd = 0.15,
    rho_sigma_prior_sd     = 0.35,
    sigma_case_prior_sd    = 0.75
  ),
  "3-component" = list(
    A0_sigma_prior_sd             = 0.50,
    gB_sigma_prior_sd             = 2.00,
    c_long_sigma_prior_sd         = 0.35,
    c_intermediate_sigma_prior_sd = 0.35,
    c_short_sigma_prior_sd        = 0.50,
    rho_sigma_prior_sd            = 0.75,
    gamma_sigma_prior_sd          = 0.35,
    sigma_case_prior_sd           = 0.75
  )
)

# Degrees of freedom of the Student-t long-rate effects (two-component model).
U_C_LONG_PRIOR_DF <- 4

# ------------------------------------------------------------------------------
# Sampler settings (as in the paper)
# ------------------------------------------------------------------------------
# QUICK_RUN shortens the chains to check that everything runs; results from a
# quick run are not reliable.

QUICK_RUN <- FALSE

SAMPLER_SETTINGS <- list(
  "2-component" = list(adapt_delta = 0.95, max_treedepth = 14),
  "3-component" = list(adapt_delta = 0.98, max_treedepth = 15)
)

N_CHAINS     <- 4
ITER_WARMUP  <- if (QUICK_RUN) 500 else 2000
ITER_SAMPLING <- if (QUICK_RUN) 500 else 4000
STAN_SEED    <- 55
STAN_INIT    <- 0.1

# ------------------------------------------------------------------------------
# Seropositivity, seroprevalence and seroreversion (03_summarise.R)
# ------------------------------------------------------------------------------

MONTHS_PER_YEAR <- 12

# Posterior draws used for the summaries: one random subset for the cohort
# summaries and another for the individual trajectories.
N_DRAWS_SUMMARY <- 200
SUMMARY_SEED    <- 2
INDIVIDUAL_SEED <- 202

# Seropositivity thresholds from the unexposed controls:
# exp(mean(log nMFI) + k * sd(log nMFI)), k = 2 and 3.
THRESHOLD_SDS <- c(
  "nMFI controls: mean+2SD" = 2,
  "nMFI controls: mean+3SD" = 3
)

# Seroreversion is measured from an antibody-specific landmark (days after
# presentation), among individuals seropositive at that landmark.
LANDMARK_DAYS <- c(IgM = 10, IgG = 14, IgG1 = 21, IgG3 = 14)

# Horizon and time grid for seroprevalence and seroreversion: steps of
# 0.1 month up to 6 months after presentation, then monthly, up to 60 years.
SERO_HORIZON_YEARS  <- 60
SERO_FINE_END_MO    <- 6
SERO_FINE_STEP_MO   <- 0.1
SERO_LATE_STEP_MO   <- 1

# Grid for the cohort median kinetics and the individual trajectories: from
# 9 months before presentation to the last sample (for individual
# trajectories, that individual's last sample), in steps of half a month.
KINETICS_START_MO <- -0.75 * MONTHS_PER_YEAR
KINETICS_STEP_MO  <- 0.5

# ------------------------------------------------------------------------------
# Figures (04_figures.R)
# ------------------------------------------------------------------------------

# Individual shown in the individual kinetics grid. NULL picks the densely
# sampled case with the most samples.
REPRESENTATIVE_ID <- NULL

ANTIBODY_COLOURS <- c(
  "IgM"  = "#3F3F3F",
  "IgG"  = "#007C89",
  "IgG1" = "#3A5DAE",
  "IgG3" = "#7A3E9D"
)

# ------------------------------------------------------------------------------
# Output folders
# ------------------------------------------------------------------------------

DIR_DATA      <- "data"
DIR_FITS      <- file.path("output", "fits")
DIR_SUMMARIES <- file.path("output", "summaries")
DIR_FIGURES   <- file.path("output", "figures")
