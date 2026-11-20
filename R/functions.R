# Functions shared by the scripts. Requires R/config.R to be sourced first.

# ------------------------------------------------------------------------------
# Antibody kinetics
# ------------------------------------------------------------------------------

# (exp(-a * dt) - exp(-b * dt)) / (b - a), with its limit dt * exp(-a * dt)
# when a and b coincide. a and b hold one value per row of the matrix dt.
diff_exp_over <- function(a, b, dt, tol = 1e-8) {
  n <- nrow(dt)
  a <- rep_len(a, n)
  b <- rep_len(b, n)
  d <- b - a
  out <- (exp(-a * dt) - exp(-b * dt)) / d
  near <- abs(d) <= tol
  if (any(near)) {
    dt_near <- dt[near, , drop = FALSE]
    out[near, ] <- dt_near * exp(-0.5 * (a[near] + b[near]) * dt_near)
  }
  out
}

# Latent antibody level for several individuals at times t_mo (months since
# presentation), as in the Stan models:
#   A(t) = A0 + gB * sum_k w_k * (exp(-c_k tau) - exp(-r tau)) / (r - c_k),
#   tau = max(t - delta, 0).
# `ind` is a list or data frame with one value per individual: A0, gB, delta,
# r, rho, c_long, c_short and, for the three-component model, gamma and
# c_intermediate. Returns a matrix with one row per individual and one column
# per time.
ab_kinetics <- function(t_mo, ind) {
  dt <- pmax(outer(-ind$delta, t_mo, "+"), 0)
  three_component <- !is.null(ind$gamma)

  w_short <- ind$rho
  w_long <- if (three_component) {
    (1 - ind$rho) * (1 - ind$gamma)
  } else {
    1 - ind$rho
  }

  k <- w_short * diff_exp_over(ind$c_short, ind$r, dt) +
    w_long * diff_exp_over(ind$c_long, ind$r, dt)

  if (three_component) {
    w_intermediate <- (1 - ind$rho) * ind$gamma
    k <- k + w_intermediate * diff_exp_over(ind$c_intermediate, ind$r, dt)
  }

  ind$A0 + ind$gB * k
}

# ------------------------------------------------------------------------------
# Simulation inputs
# ------------------------------------------------------------------------------

# Population parameters for one antigen-antibody combination, from the
# posterior medians in Supplementary Data 1, converted to the scales used in
# the Stan models (log rates and levels, logit weights).
read_population_parameters <- function(file, antigen_sel, antibody_sel) {
  s1 <- readr::read_csv(file, show_col_types = FALSE)
  s1 <- s1[s1$antigen == antigen_sel & s1$antibody == antibody_sel, ]
  if (nrow(s1) == 0) {
    stop("No parameters for ", antigen_sel, " ", antibody_sel, " in ", file)
  }
  md <- setNames(s1$md, s1$parameter)

  pars <- list(
    log_A0_pop    = log(md[["A0_pop"]]),
    log_gB_pop    = log(md[["gB_pop"]]),
    log_r         = log(md[["r_pop_mo"]]),
    log_c_long    = log(md[["c_long_pop_mo"]]),
    log_c_short   = log(md[["c_short_pop_mo"]]),
    rho_pop       = qlogis(md[["rho_pop_prob"]]),
    A0_sigma      = md[["A0_sigma"]],
    gB_sigma      = md[["gB_sigma"]],
    c_long_sigma  = md[["c_long_sigma"]],
    c_short_sigma = md[["c_short_sigma"]],
    rho_sigma     = md[["rho_sigma"]],
    sigma_case    = md[["sigma_case"]],
    sigma_ctrl    = md[["sigma_ctrl"]]
  )

  if (model_for_antibody(antibody_sel) == "3-component") {
    pars$log_c_intermediate   <- log(md[["c_intermediate_pop_mo"]])
    pars$c_intermediate_sigma <- md[["c_intermediate_sigma"]]
    pars$gamma_pop            <- qlogis(md[["gamma_pop_prob"]])
    pars$gamma_sigma          <- md[["gamma_sigma"]]
  }

  pars
}

# Kinetic parameters for n new individuals drawn from the population
# distribution. Long-rate effects are Student-t in the two-component model.
simulate_individual_parameters <- function(pars, n, model) {
  u_c_long <- if (model == "2-component") {
    rt(n, df = U_C_LONG_PRIOR_DF)
  } else {
    rnorm(n)
  }

  ind <- list(
    A0      = exp(pars$log_A0_pop + pars$A0_sigma * rnorm(n)),
    gB      = exp(pars$log_gB_pop + pars$gB_sigma * rnorm(n)),
    r       = rep(exp(pars$log_r), n),
    rho     = plogis(pars$rho_pop + pars$rho_sigma * rnorm(n)),
    c_long  = exp(pars$log_c_long + pars$c_long_sigma * u_c_long),
    c_short = exp(pars$log_c_short + pars$c_short_sigma * rnorm(n))
  )

  if (model == "3-component") {
    ind$gamma <- plogis(pars$gamma_pop + pars$gamma_sigma * rnorm(n))
    ind$c_intermediate <- exp(
      pars$log_c_intermediate + pars$c_intermediate_sigma * rnorm(n)
    )
  }

  ind
}

# ------------------------------------------------------------------------------
# Stan data and data-based priors
# ------------------------------------------------------------------------------

# Priors for the baseline and the control noise, from the controls.
make_scale_priors <- function(ctrl_y) {
  log_ctrl <- log(ctrl_y[is.finite(ctrl_y) & ctrl_y > 0])
  if (length(log_ctrl) < 3) {
    stop("Need at least 3 positive controls to set the scale priors.")
  }
  list(
    A0_prior_loc = mean(log_ctrl),
    A0_prior_sd  = max(sd(log_ctrl), 0.5),
    # Centred on an equal split of the control spread between baseline
    # variation (A0_sigma) and single-sample scatter (sigma_ctrl).
    sigma_ctrl_prior_loc = log(max(sd(log_ctrl) / sqrt(2), 0.15)),
    sigma_ctrl_prior_sd  = 0.50
  )
}

# Peak height of the kernel (exp(-c t) - exp(-r t)) / (r - c).
kernel_peak <- function(c, r) {
  if (abs(c - r) < 1e-8) {
    return(exp(-1) / c)
  }
  t_peak <- log(c / r) / (c - r)
  (exp(-c * t_peak) - exp(-r * t_peak)) / (r - c)
}

# Centre of the log_gB_pop prior ("peak rule"): the mean log early-peak excess
# of the densely sampled cases, divided by the short-component peak at the
# prior medians of rho, c_short and r. Uses the data and the priors only.
make_gB_prior_loc <- function(obs, t, nd, delta_id, A0_prior_loc,
                              kinetic_priors) {
  idx <- rep(seq_along(nd), nd)
  early <- delta_id[idx] > 0 & t <= GB_PRIOR_PEAK_WINDOW_MONTHS

  if (!any(early)) {
    warning("No early samples from densely sampled cases; ",
            "using log(50000) for the log_gB_pop prior.")
    return(log(50000))
  }

  peak <- tapply(obs[early], idx[early], max)
  baseline <- exp(A0_prior_loc)
  peak_excess <- pmax(peak - baseline, baseline)

  mean(log(peak_excess)) -
    log(
      plogis(kinetic_priors$rho_prior_loc) *
        kernel_peak(
          exp(kinetic_priors$log_c_short_prior_loc),
          exp(kinetic_priors$log_r_prior_loc)
        )
    )
}

# Stan data for one antigen-antibody combination.
# `dat` has columns sample_type ("case" or "control"), study_id, t_days,
# nmfi, is_upper_censored and upper_limit_nmfi. For censored case samples,
# nmfi holds the upper limit.
prepare_stan_data <- function(dat, antibody_sel) {
  model <- model_for_antibody(antibody_sel)
  kinetic_priors <- KINETIC_PRIORS[[antibody_sel]]

  dat <- dplyr::filter(dat, is.finite(nmfi), nmfi > 0)

  # Cases, ordered by individual and time.
  cases <- dat %>%
    dplyr::filter(sample_type == "case") %>%
    dplyr::mutate(t_months = t_days / DAYS_PER_MONTH) %>%
    dplyr::arrange(study_id, t_months)

  if (anyNA(cases$t_days) || anyNA(cases$is_upper_censored)) {
    stop("Case samples need a sampling time and a censoring indicator.")
  }
  if (any(cases$is_upper_censored == 1 &
          !(is.finite(cases$upper_limit_nmfi) & cases$upper_limit_nmfi > 0))) {
    stop("Censored case samples need a positive upper limit.")
  }
  # Stan needs a positive upper limit for every case sample; it is used only
  # for censored samples. Missing limits get the median of the others.
  has_limit <- is.finite(cases$upper_limit_nmfi) & cases$upper_limit_nmfi > 0
  if (!all(has_limit)) {
    cases$upper_limit_nmfi[!has_limit] <- if (any(has_limit)) {
      median(cases$upper_limit_nmfi[has_limit])
    } else {
      max(cases$nmfi)
    }
  }

  # Individual index; delays are estimated only for individuals with enough
  # early samples.
  case_ids <- cases %>%
    dplyr::group_by(study_id) %>%
    dplyr::summarise(
      has_delta = as.integer(
        sum(t_months <= DELTA_ESTIMATE_CUTOFF_MONTHS) >= DELTA_MIN_EARLY_OBS
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(study_id) %>%
    dplyr::mutate(
      id = dplyr::row_number(),
      delta_id = dplyr::if_else(has_delta == 1L, cumsum(has_delta), 0L)
    )

  cases <- dplyr::left_join(cases, case_ids, by = "study_id")
  nd <- as.integer(table(factor(cases$id, levels = case_ids$id)))

  controls <- dat %>%
    dplyr::filter(sample_type == "control") %>%
    dplyr::arrange(study_id)
  ctrl_ids <- dplyr::distinct(controls, study_id) %>%
    dplyr::mutate(idc = dplyr::row_number())
  controls <- dplyr::left_join(controls, ctrl_ids, by = "study_id") %>%
    dplyr::arrange(idc)
  ndc <- as.integer(table(factor(controls$idc, levels = ctrl_ids$idc)))

  # Priors.
  scale_priors <- make_scale_priors(controls$nmfi)
  log_gB_prior_loc <- make_gB_prior_loc(
    obs = cases$nmfi,
    t = cases$t_months,
    nd = nd,
    delta_id = case_ids$delta_id,
    A0_prior_loc = scale_priors$A0_prior_loc,
    kinetic_priors = kinetic_priors
  )

  stan_data <- c(
    list(
      N = nrow(cases),
      np = nrow(case_ids),
      nd = nd,
      t = cases$t_months,
      obs = cases$nmfi,
      is_upper_censored = as.integer(cases$is_upper_censored),
      upper_limit = cases$upper_limit_nmfi,
      n_delta = as.integer(max(case_ids$delta_id)),
      delta_id = as.integer(case_ids$delta_id),
      delta_lower = DELTA_RANGE[1],
      delta_upper = DELTA_RANGE[2],
      nc = nrow(controls),
      npc = nrow(ctrl_ids),
      ndc = ndc,
      obs_ctrl = controls$nmfi,
      log_gB_prior_loc = log_gB_prior_loc,
      log_gB_prior_sd = LOG_GB_PRIOR_SD,
      prior_only = 0L
    ),
    scale_priors,
    kinetic_priors,
    RE_SD_PRIORS[[model]]
  )

  if (model == "2-component") {
    stan_data$u_c_long_prior_df <- U_C_LONG_PRIOR_DF
  }

  list(
    stan_data = stan_data,
    cases = cases,
    case_ids = case_ids,
    controls = controls,
    model = model
  )
}

# ------------------------------------------------------------------------------
# Seropositivity and seroreversion
# ------------------------------------------------------------------------------

# Thresholds from the unexposed controls: exp(mean + k * SD) of log nMFI.
control_thresholds <- function(ctrl_y) {
  log_ctrl <- log(ctrl_y[is.finite(ctrl_y) & ctrl_y > 0])
  data.frame(
    threshold_label = names(THRESHOLD_SDS),
    threshold = exp(mean(log_ctrl) + THRESHOLD_SDS * sd(log_ctrl)),
    row.names = NULL
  )
}

# Time grid (months): fine steps up to fine_end_mo, monthly afterwards.
make_time_grid <- function(end_mo, fine_end_mo = SERO_FINE_END_MO) {
  fine_end_mo <- min(max(fine_end_mo, 0), end_mo)
  grid <- sort(unique(c(
    seq(0, fine_end_mo, by = SERO_FINE_STEP_MO),
    seq(fine_end_mo, end_mo, by = SERO_LATE_STEP_MO),
    end_mo
  )))
  grid[grid >= 0 & grid <= end_mo + 1e-10]
}

# Seroprevalence that never increases with time: at every time, the highest
# prevalence at that time or later.
monotone_envelope <- function(prop) {
  rev(cummax(rev(prop)))
}

# Time from the landmark to seroreversion for each individual.
# `mu` has one row per individual and one column per elapsed time `u_mo`
# (months since the landmark, starting at 0). Returns years: NA if the
# individual is not seropositive at the landmark, Inf if the level stays above
# the threshold until the end of the grid. The crossing time is interpolated
# linearly between grid points.
time_to_seroreversion <- function(u_mo, mu, threshold) {
  out <- rep(NA_real_, nrow(mu))
  positive_at_landmark <- mu[, 1] > threshold
  if (!any(positive_at_landmark)) {
    return(out)
  }

  below <- mu[, -1, drop = FALSE] <= threshold
  any_below <- rowSums(below) > 0
  out[positive_at_landmark & !any_below] <- Inf

  rows <- which(positive_at_landmark & any_below)
  if (length(rows) > 0) {
    i2 <- max.col(below[rows, , drop = FALSE], ties.method = "first") + 1L
    i1 <- i2 - 1L
    y1 <- mu[cbind(rows, i1)]
    y2 <- mu[cbind(rows, i2)]
    frac <- ifelse(y1 == y2, 1, pmin(pmax((threshold - y1) / (y2 - y1), 0), 1))
    out[rows] <- (u_mo[i1] + frac * (u_mo[i2] - u_mo[i1])) / MONTHS_PER_YEAR
  }

  out
}

# Quantile that handles infinite values (median not reached by the horizon).
quantile_type1 <- function(x, prob) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(quantile(x, probs = prob, type = 1, names = FALSE))
}

# ------------------------------------------------------------------------------
# Posterior draws
# ------------------------------------------------------------------------------

# Matrix of draws (rows) by individual (columns) for one individual-level
# variable, e.g. "A0_individual".
individual_draws <- function(draws_df, variable, np) {
  cols <- sprintf("%s[%d]", variable, seq_len(np))
  as.matrix(as.data.frame(draws_df)[, cols, drop = FALSE])
}

# Individual parameters for one posterior draw (row j), in the form used by
# ab_kinetics().
individual_parameters_for_draw <- function(ind_mats, j) {
  lapply(ind_mats, function(m) m[j, ])
}

combination_label <- function(antigen, antibody) {
  paste(antigen, antibody, sep = "_")
}
