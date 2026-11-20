# pf-ab-kinetics-longterm

Long-term antibody kinetics after *Plasmodium falciparum* malaria:
simplified analysis code and simulated data for

> Yman V, Xi et al. *[Title]*. *[Journal, year]*. doi: [DOI]

The study data cannot be shared publicly (see [Data availability](#data-availability)).
This repository instead runs the analysis on **simulated data** generated
from the posterior parameter estimates of the paper, for a cohort with the
same size and sampling design as the study (46 cases and 30 controls per
antigen–antibody combination). It runs end to end: simulate data, fit the
models and produce the main outputs.

- The two Stan models in `stan/` contain the model code used in the paper,
  unchanged. Only their header comments were shortened.
- The R code is a simplified version of the full analysis. It uses the same
  priors, sampler settings and definitions of seropositivity, seroprevalence
  and seroreversion.
- Results come from simulated data, so they will not match the numbers in
  the paper.

## Contents

| File | What it does |
|---|---|
| `01_simulate_data.R` | Simulates cases and controls for all 20 antigen–antibody combinations |
| `02_fit_models.R` | Fits the two-component (IgM) or three-component (IgG, IgG1, IgG3) model |
| `03_summarise.R` | Individual kinetics, cohort median kinetics, seroprevalence, time to seroreversion |
| `04_figures.R` | Figures from the summaries |
| `run_all.R` | Runs the four scripts in order |
| `R/config.R` | All settings: combinations, simulated design, priors, sampler settings, landmarks, thresholds |
| `R/functions.R` | Model kernel, data preparation, priors from the data, seroreversion |
| `stan/Ab_kin_2comp.stan` | Two-component model |
| `stan/Ab_kin_3comp.stan` | Three-component model |
| `data/Supplementary_Data_1_posterior_parameter_summary.csv` | Posterior estimates of the population parameters for each antigen–antibody combination (Supplementary Data 1 of the paper) |
| `data/simulated_data.csv` | Simulated data, one row per sample |

## Requirements

- R (developed with R 4.6) and the packages `dplyr`, `readr`, `ggplot2`,
  `ggh4x`, `scales`, `posterior` and `cmdstanr`.
- CmdStan, installed through cmdstanr.

```r
install.packages(c("dplyr", "readr", "ggplot2", "ggh4x", "scales", "posterior"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## How to run

Clone the repository, set the repository folder as the working directory and
run:

```r
source("run_all.R")
```

All inputs are in `data/`. `01_simulate_data.R` recreates the simulated data
exactly from the posterior estimates (fixed seed).

By default two combinations are fitted: MSP1 IgM (two-component model) and
MSP1 IgG (three-component model). In `R/config.R`:

- `FIT_ALL_COMBINATIONS <- TRUE` fits all 20 combinations;
- `COMBINATIONS_TO_FIT` selects other combinations;
- `QUICK_RUN <- TRUE` runs short chains to check that everything works.
  Results from a quick run are not reliable.

Run time with the paper's sampler settings: [to be measured].

## Simulated data

`01_simulate_data.R` simulates data for all 20 antigen–antibody combinations
(MSP1, MSP2, MSP3, AMA1 and GAMA; IgM, IgG, IgG1 and IgG3). The simulated
cohort has the size and sampling structure of the study, but not its actual
sampling times, and no study participant is reproduced.

**Parameters.** For each combination, the population parameters are the
posterior medians in `data/Supplementary_Data_1_posterior_parameter_summary.csv`:
population baseline, boost, boosting rate, decay rates, component weights,
between-individual SDs and residual SDs.

**Sampling design.** One synthetic cohort and one set of sampling times are
used for every combination:

- 14 densely sampled cases: nominal visits at 0, 7, 30, 90, 180 and 365 days
  after presentation, each shifted randomly by up to ±20%; each visit after
  day 0 is missed with probability 0.08; plus one sample 5–10 years after
  presentation.
- 32 cases sampled only long after infection: a first sample 8–45 years
  after presentation, and with probabilities 0.5, 0.35 and 0.15 one, two or
  three samples in total, 1–10 years apart and up to 45 years.
- 30 unexposed controls with one sample each.

This gives 46 cases with about 140 case samples per combination, close to the
study (about 150).

**Individuals.** For each combination, every case gets new individual
parameters drawn from the population distribution of the fitted model:

- baseline, boost and decay rates lognormal around the population values,
  with the between-individual SDs from the posterior. In the two-component
  model (IgM) the long-rate effects follow a Student-t distribution with
  4 degrees of freedom, as in the model;
- component weights (ρ, and γ in the three-component model) normal on the
  logit scale;
- the boosting rate r shared by all individuals;
- a delay δ for densely sampled cases, drawn from a normal distribution with
  SD 1/3 month truncated to ±1 month (zero for the other cases).

**Observations.** Each sample is the latent level from the model kernel
(see [Model](#model)) at the sampling time, multiplied by lognormal
measurement error with the posterior residual SD. Values at or above the upper
limit of quantification (28,000 nMFI) are set to the limit and flagged as
censored. Each control has a baseline drawn from the population baseline
distribution and one sample with the control residual SD.

**Limitations.** Levels are simulated independently for each combination, so
correlation within a person across antigens and antibodies is not
reproduced. The data cleaning, normalisation and exclusions of the study do
not apply to simulated data.

The cohort size can be changed with `N_DENSE_CASES`, `N_LATE_CASES` and
`N_CONTROLS` in `R/config.R`.

## Model

For individual *i* at time *t* (months since presentation), with delay δ and
τ = max(*t* − δ, 0):

A(*t*) = A0 + gB · Σₖ wₖ · (e^(−cₖτ) − e^(−rτ)) / (r − cₖ)

- Two-component model (IgM): a short- and a long-lived component, with
  weights ρ and 1 − ρ.
- Three-component model (IgG, IgG1, IgG3): short-, intermediate- and
  long-lived components, with weights ρ, (1 − ρ)γ and (1 − ρ)(1 − γ).

Baseline, boost, decay rates and weights vary between individuals
(hierarchical model). Observations are lognormal around A(*t*), and values at
the upper limit are treated as right-censored. Controls inform the baseline
distribution. Their individual baselines are integrated out analytically.

Seropositivity thresholds are the control mean + 2 SD and mean + 3 SD of log
nMFI. Seroprevalence and time to seroreversion are calculated from the latent
(modelled) trajectories of all modelled individuals. `03_summarise.R` gives
the exact definitions.

## Outputs

- `output/fits/`: posterior draws and convergence summaries for each fit.
- `output/summaries/`: tables behind the figures.
- `output/figures/`, in the layout of the paper's figures (PNG and PDF):
  - `individual_kinetics_grid`: one densely sampled individual in every
    fitted combination (`REPRESENTATIVE_ID` in `R/config.R`);
  - `individual_trajectories_all_individuals_<antigen>_<antibody>_<model>.pdf`:
    every modelled individual in one combination;
  - `modelled_cohort_median_and_50_95pct_individual_distribution`: median
    kinetics with the ranges covering 50% and 95% of modelled individuals;
  - `modelled_cohort_common_time_monotone_seroprevalence`: seroprevalence
    over time since presentation;
  - `modelled_cohort_seropositive_duration_antibody_landmarks`: median time
    from the landmark to seroreversion.

Kinetics are plotted on a linear time axis, seroprevalence and seroreversion
on a pseudo-log time axis. Summarising takes about 10 seconds per
combination.

## Differences from the full analysis

Not included here:

- data cleaning and normalisation of the multiplex assay data;
- exclusion of controls with high values and of one outlying IgM observation
  (described in the Methods);
- posterior predictive checks and other model diagnostics beyond R-hat,
  effective sample sizes and sampler warnings;
- comparisons between groups by prior malaria exposure;
- the two- versus three-component model comparison for IgG;
- sensitivity analyses;
- supplementary tables.

## Data availability

Due to data protection legislation and confidentiality agreements, the
pseudonymised individual-level data generated in this study cannot be made
publicly available. The data can be accessed for validation or reuse by
contacting the corresponding author. Access will be evaluated according to
the data access plan and Swedish legislation, and requires approval from the
relevant ethics and data access committees and a written data access
agreement.

Because the individual-level data cannot be shared publicly, this repository
contains the posterior parameter estimates of the paper
(`data/Supplementary_Data_1_posterior_parameter_summary.csv`), data simulated
from these estimates for a cohort with the same size and sampling design as
the study, and code to reproduce the analysis and figures on the simulated
data. Upon acceptance for publication, the repository will be archived on
Zenodo with a DOI.

## License

The code and data in this repository are released under the MIT License. See
[LICENSE](LICENSE).
