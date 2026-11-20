# 04_figures.R
#
# Figures from the summaries written by 03_summarise.R (output/summaries/),
# in the layout of the paper's figures:
#
#   individual_kinetics_grid
#       One densely sampled individual (REPRESENTATIVE_ID in R/config.R) in
#       every fitted combination: observations, posterior median latent level
#       and 95% credible interval. Dashed line: control mean+2SD threshold;
#       red dotted line: upper limit, where an observation was censored.
#   individual_trajectories_all_individuals_<antigen>_<antibody>_<model>.pdf
#       The same for every modelled individual in one combination; censored
#       observations are red triangles.
#   modelled_cohort_median_and_50_95pct_individual_distribution
#       Line: median latent level across modelled individuals. Bands: where
#       50% (25-75%) and 95% (2.5-97.5%) of the modelled individuals lie.
#       All three are posterior medians. Dashed line: control mean+2SD.
#   modelled_cohort_common_time_monotone_seroprevalence
#       Percent of modelled individuals seropositive by time since
#       presentation: posterior median (step line) and 95% CrI (band), for
#       the mean+2SD (solid) and mean+3SD (dashed) thresholds.
#   modelled_cohort_seropositive_duration_antibody_landmarks
#       Median time from the antibody-specific landmark to seroreversion
#       among individuals seropositive at the landmark: posterior median and
#       95% CrI. Times not reached within 60 years are drawn at the right edge.
#
# Kinetics are plotted on a linear time axis; seroprevalence and seroreversion
# on a pseudo-log time axis.

options(scipen = 999)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggh4x)
})

source("R/config.R")
source("R/functions.R")

dir.create(DIR_FIGURES, recursive = TRUE, showWarnings = FALSE)

read_summary <- function(name) {
  read_csv(file.path(DIR_SUMMARIES, paste0(name, ".csv")), show_col_types = FALSE) %>%
    mutate(
      antigen = factor(antigen, levels = ANTIGENS),
      antibody = factor(antibody, levels = ANTIBODIES)
    )
}

individual     <- read_summary("individual")
kinetics       <- read_summary("kinetics")
seroprevalence <- read_summary("seroprevalence")
seroreversion  <- read_summary("seroreversion")
observations   <- read_summary("observations")
thresholds     <- read_summary("thresholds")

THRESHOLD_LABELS <- names(THRESHOLD_SDS)
THRESHOLD_LINETYPES <- setNames(c("solid", "22"), THRESHOLD_LABELS)

main_threshold <- filter(thresholds, threshold_label == THRESHOLD_LABELS[1])

n_antigens   <- n_distinct(kinetics$antigen)
n_antibodies <- n_distinct(kinetics$antibody)

save_both <- function(plot, name, width_cm, height_cm, scale) {
  for (ext in c("png", "pdf")) {
    ggsave(
      file.path(DIR_FIGURES, paste0(name, ".", ext)), plot,
      width = width_cm, height = height_cm, units = "cm",
      scale = scale, dpi = 300
    )
  }
}

# Upper-limit lines, only where an observation was censored.
upper_limit_lines <- function(obs, ...) {
  obs %>%
    filter(is_upper_censored == 1, is.finite(upper_limit_nmfi)) %>%
    distinct(..., upper_limit = upper_limit_nmfi)
}

# ------------------------------------------------------------------------------
# 1. Individual kinetics: one individual across all fitted combinations
# ------------------------------------------------------------------------------

representative_id <- REPRESENTATIVE_ID
if (is.null(representative_id)) {
  representative_id <- observations %>%
    filter(has_delta == 1) %>%
    count(study_id) %>%
    arrange(desc(n), study_id) %>%
    slice(1) %>%
    pull(study_id)
}

rep_pred <- filter(individual, study_id == representative_id)
rep_obs  <- filter(observations, study_id == representative_id)

p_ind <- ggplot() +
  geom_hline(
    data = upper_limit_lines(rep_obs, antigen, antibody),
    aes(yintercept = upper_limit),
    linewidth = 0.35, linetype = "dotted", colour = "red", alpha = 0.55
  ) +
  geom_hline(
    data = main_threshold,
    aes(yintercept = threshold),
    linewidth = 0.55, linetype = "dashed", colour = "black"
  ) +
  geom_ribbon(
    data = rep_pred,
    aes(x = t_years, ymin = lo, ymax = hi, fill = antibody),
    alpha = 0.1, show.legend = FALSE
  ) +
  geom_line(
    data = rep_pred,
    aes(x = t_years, y = md, colour = antibody),
    linewidth = 0.7
  ) +
  geom_point(
    data = rep_obs,
    aes(x = t_years, y = nmfi),
    size = 2, colour = "black"
  ) +
  scale_y_log10() +
  scale_fill_manual(values = ANTIBODY_COLOURS) +
  scale_colour_manual(values = ANTIBODY_COLOURS) +
  labs(x = "Years since presentation", y = "Antibody level (nMFI)") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    panel.spacing.y = grid::unit(0.5, "cm")
  ) +
  facet_grid2(
    rows = vars(antibody), cols = vars(antigen),
    scales = "free_y", independent = "y", axes = "y", remove_labels = "x"
  )

save_both(
  p_ind, "individual_kinetics_grid",
  width_cm = 4 * n_antigens, height_cm = 1.5 + 3.375 * n_antibodies, scale = 2
)

# ------------------------------------------------------------------------------
# 2. Individual kinetics: all modelled individuals, one file per combination
# ------------------------------------------------------------------------------

for (combo in split(distinct(individual, antigen, antibody, model),
                    seq_len(nrow(distinct(individual, antigen, antibody))))) {
  d_tr <- semi_join(individual, combo, by = c("antigen", "antibody"))
  d_ob <- semi_join(observations, combo, by = c("antigen", "antibody"))

  id_levels <- unique(d_tr$study_id)
  d_tr$study_id <- factor(d_tr$study_id, levels = id_levels)
  d_ob$study_id <- factor(d_ob$study_id, levels = id_levels)
  n_rows <- ceiling(length(id_levels) / 5)

  p <- ggplot() +
    geom_hline(
      data = upper_limit_lines(d_ob, study_id),
      aes(yintercept = upper_limit),
      linewidth = 0.30, linetype = "dotted", colour = "red", alpha = 0.55
    ) +
    geom_ribbon(data = d_tr, aes(x = t_years, ymin = lo, ymax = hi), alpha = 0.10) +
    geom_line(data = d_tr, aes(x = t_years, y = md), linewidth = 0.6) +
    geom_point(
      data = filter(d_ob, is_upper_censored == 0),
      aes(x = t_years, y = nmfi), size = 0.9
    ) +
    geom_point(
      data = filter(d_ob, is_upper_censored == 1),
      aes(x = t_years, y = nmfi),
      size = 1.2, shape = 24, fill = "red", colour = "black"
    ) +
    scale_y_log10() +
    facet_wrap(~ study_id, scales = "free", ncol = 5) +
    labs(
      x = "Years since presentation",
      y = "nMFI antibody level, log scale",
      title = sprintf(
        "%s %s, individual trajectories, all modelled individuals",
        as.character(combo$antigen), as.character(combo$antibody)
      )
    ) +
    theme_bw(base_size = 9) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 10)
    )

  ggsave(
    file.path(DIR_FIGURES, sprintf(
      "individual_trajectories_all_individuals_%s_%s_%s.pdf",
      as.character(combo$antigen), as.character(combo$antibody), combo$model
    )),
    p, width = 10, height = max(8, 1.1 * n_rows), scale = 1.2
  )
}

# ------------------------------------------------------------------------------
# 3. Cohort median kinetics with 50% and 95% of modelled individuals
# ------------------------------------------------------------------------------

p_cohort <- ggplot(kinetics, aes(x = t_years)) +
  geom_ribbon(
    aes(ymin = q2.5, ymax = q97.5, fill = antibody),
    alpha = 0.08, colour = NA, show.legend = FALSE
  ) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75, fill = antibody),
    alpha = 0.22, colour = NA, show.legend = FALSE
  ) +
  geom_hline(
    data = main_threshold,
    aes(yintercept = threshold),
    linewidth = 0.55, colour = "black", linetype = "dashed"
  ) +
  geom_line(aes(y = q50, colour = antibody), linewidth = 0.8) +
  scale_fill_manual(values = ANTIBODY_COLOURS) +
  scale_colour_manual(values = ANTIBODY_COLOURS) +
  scale_y_log10() +
  labs(x = "Years since presentation", y = "Antibody level (nMFI)") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  facet_grid2(
    rows = vars(antibody), cols = vars(antigen),
    scales = "free_y", independent = "y"
  )

save_both(
  p_cohort, "modelled_cohort_median_and_50_95pct_individual_distribution",
  width_cm = 4 * n_antigens, height_cm = 1.5 + 3.375 * n_antibodies, scale = 2
)

# ------------------------------------------------------------------------------
# 4. Seroprevalence over time since presentation
# ------------------------------------------------------------------------------

SERO_TIME_BREAKS <- c(0, 2 / 12, 0.5, 1, 2, 5, 10, 20, 40)
SERO_TIME_LABELS <- c("0", "2mo", "6mo", "1y", "2y", "5y", "10y", "20y", "40y")
SERO_PLOT_MAX_YEARS <- 60
DURATION_PLOT_MAX_YEARS <- 70

sero_time_scale <- function(expand_mult = c(0.01, 0.03)) {
  scale_x_continuous(
    transform = scales::transform_pseudo_log(sigma = 0.1, base = 10),
    breaks = SERO_TIME_BREAKS,
    labels = SERO_TIME_LABELS,
    expand = expansion(mult = expand_mult)
  )
}

sero_plot_data <- seroprevalence %>%
  filter(t_years >= 0, t_years <= SERO_PLOT_MAX_YEARS) %>%
  mutate(threshold_label = factor(threshold_label, levels = THRESHOLD_LABELS))

p_sero <- ggplot(
  sero_plot_data,
  aes(x = t_years, group = interaction(antibody, threshold_label))
) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi, fill = antibody),
    alpha = 0.08, colour = NA, show.legend = FALSE
  ) +
  geom_step(
    aes(y = md, colour = antibody, linetype = threshold_label),
    direction = "hv", linewidth = 0.7
  ) +
  scale_colour_manual(values = ANTIBODY_COLOURS) +
  scale_fill_manual(values = ANTIBODY_COLOURS) +
  scale_linetype_manual(values = THRESHOLD_LINETYPES) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  sero_time_scale() +
  coord_cartesian(xlim = c(0, SERO_PLOT_MAX_YEARS)) +
  guides(
    fill = "none",
    colour = guide_legend(
      title = NULL,
      override.aes = list(linewidth = 1.2, linetype = "solid")
    ),
    linetype = guide_legend(
      title = "Threshold",
      override.aes = list(colour = "black", linewidth = 1.0)
    )
  ) +
  labs(x = "Time since presentation", y = "Percent seropositive") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.text.x = element_text(size = 9),
    panel.spacing.y = grid::unit(0.5, "cm"),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14)
  ) +
  facet_grid2(
    rows = vars(antibody), cols = vars(antigen),
    scales = "fixed", axes = "y", remove_labels = "x"
  )

save_both(
  p_sero, "modelled_cohort_common_time_monotone_seroprevalence",
  width_cm = 4 * n_antigens, height_cm = 3.75 * n_antibodies, scale = 1.9
)

# ------------------------------------------------------------------------------
# 5. Median time to seroreversion from the antibody-specific landmark
# ------------------------------------------------------------------------------

duration_plot_data <- seroreversion %>%
  mutate(threshold_label = factor(threshold_label, levels = THRESHOLD_LABELS)) %>%
  filter(!is.na(median_years_md))

p_duration <- ggplot(duration_plot_data, aes(y = antigen, colour = antibody)) +
  geom_linerange(aes(xmin = median_years_lo, xmax = median_years_hi), linewidth = 0.55) +
  geom_point(aes(x = median_years_md), size = 4.4) +
  scale_x_continuous(
    transform = scales::transform_pseudo_log(sigma = 0.1, base = 10),
    breaks = SERO_TIME_BREAKS,
    labels = SERO_TIME_LABELS,
    expand = expansion(mult = 0.01)
  ) +
  coord_cartesian(xlim = c(0, DURATION_PLOT_MAX_YEARS)) +
  scale_y_discrete(limits = rev) +
  scale_colour_manual(values = ANTIBODY_COLOURS) +
  labs(x = "Median time to seroreversion", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  facet_grid(rows = vars(threshold_label), cols = vars(antibody))

save_both(
  p_duration, "modelled_cohort_seropositive_duration_antibody_landmarks",
  width_cm = 4 + 3.75 * n_antibodies, height_cm = 3.5 + 1.2 * n_antigens, scale = 2
)

message("Figures written to ", DIR_FIGURES)
