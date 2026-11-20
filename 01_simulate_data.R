# 01_simulate_data.R
#
# Simulates a synthetic data set with the structure of the study, for all 20
# antigen-antibody combinations.
#
# - Sampling design: densely sampled cases (first year after presentation plus
#   one sample 5-10 years later), cases sampled only 8-45 years after
#   infection, and unexposed controls with one sample each (see R/config.R).
#   The same synthetic individuals and sampling times are used for every
#   combination.
# - Antibody levels: new individuals are drawn from the population
#   distribution of the fitted models, using the posterior medians of the
#   population parameters in the paper's Supplementary Data 1. No individual
#   from the study is reproduced.
# - Levels within an individual are simulated independently for each
#   combination, so correlation across antigens and antibodies is not
#   reproduced.
#
# Output:
#   data/simulated_data.csv   one row per sample

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/config.R")
source("R/functions.R")

if (!file.exists(SIM_PARAMETER_FILE)) {
  stop(
    "Population parameters not found: ", SIM_PARAMETER_FILE, "\n",
    "Download Supplementary Data 1 from the paper and save it there."
  )
}

# ------------------------------------------------------------------------------
# Sampling design (shared by all combinations)
# ------------------------------------------------------------------------------

set.seed(SIM_SEED)

dense_ids   <- sprintf("D%02d", seq_len(N_DENSE_CASES))
late_ids    <- sprintf("L%02d", seq_len(N_LATE_CASES))
control_ids <- sprintf("C%02d", seq_len(N_CONTROLS))

design_dense <- bind_rows(lapply(dense_ids, function(id) {
  attended <- c(TRUE, runif(length(DENSE_VISIT_DAYS) - 1) > P_MISSED_VISIT)
  days <- DENSE_VISIT_DAYS[attended]
  days <- round(days * runif(
    length(days), 1 - DENSE_VISIT_JITTER, 1 + DENSE_VISIT_JITTER
  ))
  late_day <- round(
    runif(1, DENSE_LATE_YEARS[1], DENSE_LATE_YEARS[2]) * DAYS_PER_YEAR
  )
  data.frame(study_id = id, t_days = c(days, late_day))
}))

design_late <- bind_rows(lapply(late_ids, function(id) {
  n <- sample(seq_along(LATE_N_SAMPLES_PROBS), 1, prob = LATE_N_SAMPLES_PROBS)
  first <- runif(1, LATE_FIRST_YEARS[1], LATE_FIRST_YEARS[2])
  gaps <- runif(n - 1, LATE_GAP_YEARS[1], LATE_GAP_YEARS[2])
  years <- first + c(0, cumsum(gaps))
  years <- years[years <= MAX_FOLLOW_UP_YEARS]
  data.frame(study_id = id, t_days = round(years * DAYS_PER_YEAR))
}))

design_cases <- bind_rows(design_dense, design_late)
case_ids <- c(dense_ids, late_ids)

# Delays (months): drawn for densely sampled cases from a normal distribution
# with SD DELTA_SD truncated to DELTA_RANGE; zero for the others.
p_range <- pnorm(DELTA_RANGE / DELTA_SD)
delta_dense <- qnorm(runif(N_DENSE_CASES, p_range[1], p_range[2])) * DELTA_SD
true_delta <- setNames(c(delta_dense, rep(0, N_LATE_CASES)), case_ids)

# ------------------------------------------------------------------------------
# Antibody levels for each combination
# ------------------------------------------------------------------------------

all_combinations <- expand.grid(
  antigen = ANTIGENS,
  antibody = ANTIBODIES,
  stringsAsFactors = FALSE
)

sim_data <- list()

for (k in seq_len(nrow(all_combinations))) {
  antigen_sel  <- all_combinations$antigen[k]
  antibody_sel <- all_combinations$antibody[k]
  model <- model_for_antibody(antibody_sel)

  # Separate seed per combination, so that changing one does not change the
  # others.
  set.seed(SIM_SEED + k)

  pars <- read_population_parameters(
    SIM_PARAMETER_FILE, antigen_sel, antibody_sel
  )

  # Cases.
  ind <- simulate_individual_parameters(pars, length(case_ids), model)
  ind$delta <- unname(true_delta)

  case_rows <- bind_rows(lapply(seq_along(case_ids), function(i) {
    t_days <- design_cases$t_days[design_cases$study_id == case_ids[i]]
    ind_i <- lapply(ind, `[`, i)
    mu <- as.numeric(ab_kinetics(t_days / DAYS_PER_MONTH, ind_i))
    y <- mu * exp(rnorm(length(mu), 0, pars$sigma_case))
    data.frame(
      sample_type = "case",
      study_id = case_ids[i],
      t_days = t_days,
      latent_nmfi = mu,
      nmfi = pmin(y, UPPER_LIMIT_NMFI),
      is_upper_censored = as.integer(y >= UPPER_LIMIT_NMFI),
      upper_limit_nmfi = UPPER_LIMIT_NMFI
    )
  }))

  # Controls: one sample each around an individual baseline.
  log_A0_ctrl <- pars$log_A0_pop + pars$A0_sigma * rnorm(N_CONTROLS)
  y_ctrl <- exp(log_A0_ctrl + pars$sigma_ctrl * rnorm(N_CONTROLS))
  control_rows <- data.frame(
    sample_type = "control",
    study_id = control_ids,
    t_days = NA_real_,
    latent_nmfi = exp(log_A0_ctrl),
    nmfi = y_ctrl,
    is_upper_censored = 0L,
    upper_limit_nmfi = UPPER_LIMIT_NMFI
  )

  sim_data[[k]] <- bind_rows(case_rows, control_rows) %>%
    mutate(antigen = antigen_sel, antibody = antibody_sel, .before = 1)
}

sim_data <- bind_rows(sim_data)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

dir.create(DIR_DATA, showWarnings = FALSE)

write_csv(sim_data, file.path(DIR_DATA, "simulated_data.csv"))

n_cases <- n_distinct(sim_data$study_id[sim_data$sample_type == "case"])
n_samples <- nrow(design_cases)
message(
  "Simulated ", nrow(all_combinations), " combinations: ",
  n_cases, " cases with ", n_samples, " samples and ",
  N_CONTROLS, " controls each.\n",
  "Censored at the upper limit: ",
  sum(sim_data$is_upper_censored), " of ",
  sum(sim_data$sample_type == "case"), " case samples."
)
