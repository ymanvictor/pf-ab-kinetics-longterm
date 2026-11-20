# 03_summarise.R
#
# Summarises each fit in output/fits/ for the modelled individuals (all cases
# in the fit):
#
#   1. Individual kinetics: posterior median and 95% credible interval of each
#      individual's latent antibody level, from 9 months before presentation
#      to their last sample.
#   2. Cohort median kinetics: at each time, the median and the 25-75% and
#      2.5-97.5% ranges of the latent levels across individuals; the curves
#      shown are posterior medians of these quantities.
#   3. Seroprevalence over time: the proportion of individuals whose latent
#      level is above the control threshold at each time since presentation.
#      Within each posterior draw the curve is made non-increasing (at each
#      time, the highest proportion at that time or later).
#   4. Time to seroreversion: among individuals seropositive at the
#      antibody-specific landmark (IgM 10 d, IgG 14 d, IgG1 21 d, IgG3 14 d),
#      the median time from the landmark until the latent level first falls
#      to the threshold. Individuals still above the threshold at 60 years
#      count as not reverted (infinite time).
#
# Thresholds are exp(mean + 2 SD) and exp(mean + 3 SD) of log nMFI in the
# unexposed controls. All summaries use N_DRAWS_SUMMARY posterior draws.
#
# Output: CSV files in output/summaries/.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(posterior)
})

source("R/config.R")
source("R/functions.R")

dir.create(DIR_SUMMARIES, recursive = TRUE, showWarnings = FALSE)

# Individual-level variables in the fits, by the names used in ab_kinetics().
INDIVIDUAL_VARIABLE_NAMES <- c(
  A0             = "A0_individual",
  gB             = "gB_individual",
  delta          = "delta",
  r              = "r_individual",
  rho            = "rho_individual",
  gamma          = "gamma_individual",
  c_long         = "c_long_individual",
  c_intermediate = "c_intermediate_individual",
  c_short        = "c_short_individual"
)

KINETICS_PROBS <- c(0.025, 0.25, 0.5, 0.75, 0.975)

summarise_fit <- function(fit_file) {
  f <- readRDS(fit_file)
  label <- combination_label(f$antigen, f$antibody)
  message("Summarising ", label)

  keys <- data.frame(antigen = f$antigen, antibody = f$antibody, model = f$model)
  np <- nrow(f$case_ids)
  thresholds <- control_thresholds(f$controls$nmfi)

  in_fit <- unique(sub("\\[.*$", "", variables(f$draws_individual)))
  vars <- INDIVIDUAL_VARIABLE_NAMES[INDIVIDUAL_VARIABLE_NAMES %in% in_fit]
  ind_mats_all <- lapply(vars, function(v) {
    individual_draws(f$draws_individual, v, np)
  })
  n_draws <- nrow(f$draws_individual)

  # 1. Individual kinetics: each individual from 9 months before presentation
  # to their last sample, over a random subset of draws.
  set.seed(INDIVIDUAL_SEED)
  draw_idx_ind <- sort(sample(seq_len(n_draws), min(N_DRAWS_SUMMARY, n_draws)))

  individual <- bind_rows(lapply(seq_len(np), function(i) {
    t_max_i <- max(f$cases$t_months[f$cases$id == i])
    grid_i <- seq(KINETICS_START_MO, t_max_i, by = KINETICS_STEP_MO)
    # One row per posterior draw of individual i.
    pars_i <- lapply(ind_mats_all, function(m) m[draw_idx_ind, i])
    mu <- ab_kinetics(grid_i, pars_i)
    data.frame(
      study_id = f$case_ids$study_id[i],
      has_delta = f$case_ids$has_delta[i],
      t_years = grid_i / MONTHS_PER_YEAR,
      lo = apply(mu, 2, quantile, 0.025),
      md = apply(mu, 2, median),
      hi = apply(mu, 2, quantile, 0.975)
    )
  }))

  # Posterior draws used for the cohort summaries (2-4).
  set.seed(SUMMARY_SEED)
  draw_idx <- sort(sample(seq_len(n_draws), min(N_DRAWS_SUMMARY, n_draws)))
  n_use <- length(draw_idx)
  ind_mats <- lapply(ind_mats_all, function(m) m[draw_idx, , drop = FALSE])

  # Time grids (months since presentation).
  t_max_mo <- max(f$cases$t_months)
  grid_kinetics <- seq(KINETICS_START_MO, t_max_mo, by = KINETICS_STEP_MO)
  grid_sero <- make_time_grid(SERO_HORIZON_YEARS * MONTHS_PER_YEAR)

  landmark_mo <- LANDMARK_DAYS[[f$antibody]] / DAYS_PER_MONTH
  grid_since_landmark <- make_time_grid(
    SERO_HORIZON_YEARS * MONTHS_PER_YEAR,
    fine_end_mo = max(SERO_FINE_END_MO - landmark_mo, 0)
  )
  grid_duration <- landmark_mo + grid_since_landmark

  n_thr <- nrow(thresholds)
  kinetics_q <- array(NA_real_, c(n_use, length(KINETICS_PROBS), length(grid_kinetics)))
  sero_prop <- array(NA_real_, c(n_use, n_thr, length(grid_sero)))
  reversion <- list()

  for (j in seq_len(n_use)) {
    ind <- individual_parameters_for_draw(ind_mats, j)

    # Cohort distribution of latent levels.
    mu <- ab_kinetics(grid_kinetics, ind)
    kinetics_q[j, , ] <- apply(mu, 2, quantile, probs = KINETICS_PROBS, names = FALSE)

    # Seroprevalence and time to seroreversion.
    mu_sero <- ab_kinetics(grid_sero, ind)
    mu_duration <- ab_kinetics(grid_duration, ind)

    for (h in seq_len(n_thr)) {
      thr <- thresholds$threshold[h]
      sero_prop[j, h, ] <- monotone_envelope(colMeans(mu_sero > thr))

      duration <- time_to_seroreversion(grid_since_landmark, mu_duration, thr)
      positive <- !is.na(duration)
      reversion[[length(reversion) + 1]] <- data.frame(
        threshold_label = thresholds$threshold_label[h],
        draw = j,
        prop_seropositive_at_landmark = mean(positive),
        n_seropositive_at_landmark = sum(positive),
        median_years = if (any(positive)) {
          quantile_type1(duration[positive], 0.5)
        } else {
          NA_real_
        }
      )
    }
  }

  # 2. Cohort median kinetics.
  kinetics_md <- apply(kinetics_q, c(2, 3), median)
  kinetics <- data.frame(
    t_years = grid_kinetics / MONTHS_PER_YEAR,
    q2.5  = kinetics_md[1, ],
    q25   = kinetics_md[2, ],
    q50   = kinetics_md[3, ],
    q75   = kinetics_md[4, ],
    q97.5 = kinetics_md[5, ],
    q50_lo = apply(kinetics_q[, 3, , drop = FALSE], 3, quantile, 0.025),
    q50_hi = apply(kinetics_q[, 3, , drop = FALSE], 3, quantile, 0.975)
  )

  # 3. Seroprevalence over time.
  seroprevalence <- bind_rows(lapply(seq_len(n_thr), function(h) {
    p <- sero_prop[, h, , drop = FALSE]
    data.frame(
      threshold_label = thresholds$threshold_label[h],
      t_years = grid_sero / MONTHS_PER_YEAR,
      md = apply(p, 3, median),
      lo = apply(p, 3, quantile, 0.025),
      hi = apply(p, 3, quantile, 0.975)
    )
  }))

  # 4. Time to seroreversion.
  seroreversion <- bind_rows(reversion) %>%
    group_by(threshold_label) %>%
    summarise(
      landmark_days = LANDMARK_DAYS[[f$antibody]],
      seropositive_at_landmark_md = median(prop_seropositive_at_landmark),
      seropositive_at_landmark_lo = quantile_type1(prop_seropositive_at_landmark, 0.025),
      seropositive_at_landmark_hi = quantile_type1(prop_seropositive_at_landmark, 0.975),
      median_years_md = quantile_type1(median_years[n_seropositive_at_landmark > 0], 0.5),
      median_years_lo = quantile_type1(median_years[n_seropositive_at_landmark > 0], 0.025),
      median_years_hi = quantile_type1(median_years[n_seropositive_at_landmark > 0], 0.975),
      prop_draws_median_reached = mean(is.finite(median_years[n_seropositive_at_landmark > 0])),
      horizon_years = SERO_HORIZON_YEARS,
      .groups = "drop"
    )

  # Observations as fitted, for the figures.
  observations <- f$cases %>%
    transmute(
      study_id, has_delta,
      t_years = t_months / MONTHS_PER_YEAR,
      nmfi, is_upper_censored, upper_limit_nmfi
    )

  list(
    individual     = cbind(keys, individual),
    kinetics       = cbind(keys, kinetics),
    seroprevalence = cbind(keys, seroprevalence),
    seroreversion  = cbind(keys, seroreversion),
    observations   = cbind(keys, observations),
    thresholds     = cbind(keys, thresholds)
  )
}

fit_files <- list.files(DIR_FITS, pattern = "^fit_.*\\.rds$", full.names = TRUE)
if (length(fit_files) == 0) {
  stop("No fits found in ", DIR_FITS, ". Run 02_fit_models.R first.")
}

summaries <- lapply(fit_files, summarise_fit)

# Combine across fits and save.
for (name in names(summaries[[1]])) {
  combined <- bind_rows(lapply(summaries, `[[`, name))
  write_csv(combined, file.path(DIR_SUMMARIES, paste0(name, ".csv")))
}

message("Summaries written to ", DIR_SUMMARIES)
