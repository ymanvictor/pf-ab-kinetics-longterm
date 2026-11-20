# 02_fit_models.R
#
# Fits the two-component model (IgM) or the three-component model (IgG, IgG1,
# IgG3) to each combination in COMBINATIONS_TO_FIT (R/config.R), with the
# priors and sampler settings used in the paper.
#
# Priors partly depend on the data, as in the paper:
#   - baseline and control-noise priors from the controls;
#   - the log_gB_pop prior from the early peaks of the densely sampled cases
#     ("peak rule"; see make_gB_prior_loc() in R/functions.R).
#
# Output, per combination (output/fits/):
#   fit_<antigen>_<antibody>.rds    posterior draws and the data as fitted
#   convergence_<antigen>_<antibody>.csv
#
# Run time with the paper's settings (4 chains x 6,000 iterations): minutes
# for the two-component model, longer for the three-component model. Set
# QUICK_RUN <- TRUE in R/config.R for a fast check that everything runs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(cmdstanr)
  library(posterior)
})

source("R/config.R")
source("R/functions.R")

dir.create(DIR_FITS, recursive = TRUE, showWarnings = FALSE)

sim_data <- read_csv(
  file.path(DIR_DATA, "simulated_data.csv"),
  show_col_types = FALSE
)

# Variables kept from each fit.
POPULATION_VARIABLES <- list(
  "2-component" = c(
    "log_A0_pop", "log_gB_pop", "log_r", "log_c_pop", "rho_pop",
    "A0_sigma", "gB_sigma", "c_long_sigma", "c_short_sigma", "rho_sigma",
    "sigma_case", "sigma_ctrl"
  ),
  "3-component" = c(
    "log_A0_pop", "log_gB_pop", "log_r", "log_c_pop", "rho_pop", "gamma_pop",
    "A0_sigma", "gB_sigma", "c_long_sigma", "c_intermediate_sigma",
    "c_short_sigma", "rho_sigma", "gamma_sigma", "sigma_case", "sigma_ctrl"
  )
)

INDIVIDUAL_VARIABLES <- list(
  "2-component" = c(
    "A0_individual", "gB_individual", "delta", "r_individual",
    "rho_individual", "c_long_individual", "c_short_individual"
  ),
  "3-component" = c(
    "A0_individual", "gB_individual", "delta", "r_individual",
    "rho_individual", "gamma_individual", "c_long_individual",
    "c_intermediate_individual", "c_short_individual"
  )
)

# Compile only the models that are needed.
models_needed <- unique(model_for_antibody(COMBINATIONS_TO_FIT$antibody))
stan_models <- lapply(
  setNames(models_needed, models_needed),
  function(m) cmdstan_model(STAN_FILES[[m]])
)

for (k in seq_len(nrow(COMBINATIONS_TO_FIT))) {
  antigen_sel  <- COMBINATIONS_TO_FIT$antigen[k]
  antibody_sel <- COMBINATIONS_TO_FIT$antibody[k]
  label <- combination_label(antigen_sel, antibody_sel)
  message("\n=== Fitting ", label, " ===")

  dat <- filter(sim_data, antigen == antigen_sel, antibody == antibody_sel)
  prep <- prepare_stan_data(dat, antibody_sel)
  model <- prep$model

  message(
    "  ", model, " model: ", prep$stan_data$np, " cases (",
    prep$stan_data$n_delta, " with an estimated delay), ",
    prep$stan_data$N, " case samples, ", prep$stan_data$npc, " controls"
  )

  fit <- stan_models[[model]]$sample(
    data = prep$stan_data,
    chains = N_CHAINS,
    parallel_chains = N_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    seed = STAN_SEED,
    refresh = 1000,
    adapt_delta = SAMPLER_SETTINGS[[model]]$adapt_delta,
    max_treedepth = SAMPLER_SETTINGS[[model]]$max_treedepth,
    init = STAN_INIT
  )

  # Convergence: R-hat and effective sample sizes of the population-level
  # parameters, and sampler diagnostics.
  convergence <- fit$summary(
    variables = POPULATION_VARIABLES[[model]],
    "median", "quantile2", "rhat", "ess_bulk", "ess_tail"
  )
  sampler <- fit$diagnostic_summary(quiet = TRUE)

  write_csv(
    convergence,
    file.path(DIR_FITS, paste0("convergence_", label, ".csv"))
  )

  message(
    "  Divergences: ", sum(sampler$num_divergent),
    "; max tree depth hits: ", sum(sampler$num_max_treedepth),
    "; min E-BFMI: ", round(min(sampler$ebfmi), 2),
    "; max R-hat: ", round(max(convergence$rhat, na.rm = TRUE), 3),
    "; min bulk ESS: ", round(min(convergence$ess_bulk, na.rm = TRUE))
  )

  saveRDS(
    list(
      antigen = antigen_sel,
      antibody = antibody_sel,
      model = model,
      draws_population = as_draws_df(
        fit$draws(variables = POPULATION_VARIABLES[[model]])
      ),
      draws_individual = as_draws_df(
        fit$draws(variables = INDIVIDUAL_VARIABLES[[model]])
      ),
      cases = prep$cases,
      case_ids = prep$case_ids,
      controls = prep$controls,
      stan_data = prep$stan_data,
      sampler = sampler
    ),
    file.path(DIR_FITS, paste0("fit_", label, ".rds"))
  )

  rm(fit)
  invisible(gc())
}
