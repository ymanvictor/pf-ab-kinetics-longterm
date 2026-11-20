// Two-component antibody kinetics model (used for IgM).
//
// Model code identical to the version used in the paper; only this header
// comment has been shortened.
//
// Kernel, with tau = max(t - delta, 0) in months:
//   A(t) = A0 + gB * [ w_short * k(c_short) + w_long * k(c_long) ]
//   k(c) = (exp(-c * tau) - exp(-r * tau)) / (r - c)
//   w_short = rho, w_long = 1 - rho
//
// Individual effects are normal and non-centred, except the long-rate effects
// (u_c_long), which follow a Student-t distribution with u_c_long_prior_df
// degrees of freedom (4 in the paper), so that a few individuals with
// persistent IgM can sit in the tail without inflating c_long_sigma.
//
// Control baselines are integrated out analytically; upper-censored case
// observations use the right-censored lognormal likelihood.
//
// Internal ordering: log_c_pop[1] = log c_long, log_c_pop[2] = log c_short.

functions {

  // (exp(-a * dt) - exp(-b * dt)) / (b - a), stable for b close to a.
  // The smaller rate is factored out so expm1 never overflows when a > b.
  real diff_exp_over(real a, real b, real dt) {
    real d = abs(b - a);

    if (d > 1e-12) {
      return -exp(-fmin(a, b) * dt) * expm1(-d * dt) / d;
    } else {
      return dt * exp(-0.5 * (a + b) * dt);
    }
  }

  vector AbKin(vector t, real A0, real gB, real delta, real rho,
               real c_long, real r, real c_short) {
    int T = size(t);
    vector[T] y;

    for (i in 1:T) {
      real dt = fmax(t[i] - delta, 0.0);

      y[i] = A0
        + gB * (  rho       * diff_exp_over(c_short, r, dt)
                + (1 - rho) * diff_exp_over(c_long, r, dt) );
    }

    return y;
  }
}

data {
  // --------------------------------------------------------------------------
  // CASES
  // --------------------------------------------------------------------------

  int<lower=0> N;
  int<lower=0> np;
  array[np] int<lower=0> nd;

  vector<lower=0>[N] t;

  vector<lower=0>[N] obs;

  array[N] int<lower=0, upper=1> is_upper_censored;
  vector<lower=0>[N] upper_limit;

  // Delay subset.
  int<lower=0> n_delta;
  array[np] int<lower=0, upper=n_delta> delta_id;
  real delta_lower;
  real delta_upper;

  // --------------------------------------------------------------------------
  // CONTROLS
  // --------------------------------------------------------------------------

  int<lower=0> nc;
  int<lower=0> npc;
  array[npc] int<lower=0> ndc;
  vector<lower=0>[nc] obs_ctrl;

  // --------------------------------------------------------------------------
  // Data-adaptive priors on log response / log noise scale
  // --------------------------------------------------------------------------

  real A0_prior_loc;
  real<lower=0> A0_prior_sd;

  real sigma_ctrl_prior_loc;
  real<lower=0> sigma_ctrl_prior_sd;

  real log_gB_prior_loc;
  real<lower=0> log_gB_prior_sd;

  // --------------------------------------------------------------------------
  // Kinetic priors supplied from R
  // --------------------------------------------------------------------------

  real log_r_prior_loc;
  real<lower=0> log_r_prior_sd;

  real log_c_long_prior_loc;
  real<lower=0> log_c_long_prior_sd;

  real log_c_short_prior_loc;
  real<lower=0> log_c_short_prior_sd;

  real rho_prior_loc;
  real<lower=0> rho_prior_sd;

  // --------------------------------------------------------------------------
  // Half-normal scales for random-effect SDs and case noise
  // --------------------------------------------------------------------------

  real<lower=0> A0_sigma_prior_sd;
  real<lower=0> gB_sigma_prior_sd;
  real<lower=0> c_long_sigma_prior_sd;
  real<lower=0> c_short_sigma_prior_sd;
  real<lower=0> rho_sigma_prior_sd;
  real<lower=0> sigma_case_prior_sd;

  // Degrees of freedom of the Student-t distribution of u_c_long.
  real<lower=1> u_c_long_prior_df;

  int<lower=0, upper=1> prior_only;
}

transformed data {
  // Case observation index sets for the vectorised likelihood.
  int n_obs = N - sum(is_upper_censored);
  int n_cens = sum(is_upper_censored);
  array[n_obs] int idx_obs;
  array[n_cens] int idx_cens;

  // Control sufficient statistics on the log scale.
  vector[npc] ctrl_n;
  vector[npc] ctrl_logy_mean;
  vector[npc] ctrl_logy_ss;
  // log_lik constant per control: -(n-1)/2 log(2 pi) - 0.5 log(n) - sum(log y)
  vector[npc] ctrl_ll_const;

  {
    int io = 1;
    int ic = 1;

    for (n in 1:N) {
      if (is_upper_censored[n] == 1) {
        idx_cens[ic] = n;
        ic += 1;
      } else {
        idx_obs[io] = n;
        io += 1;
      }
    }
  }

  {
    int nddc = 1;

    for (j in 1:npc) {
      int steps = ndc[j];
      vector[steps] ly = log(obs_ctrl[nddc:(nddc + steps - 1)]);

      ctrl_n[j] = steps;
      ctrl_logy_mean[j] = mean(ly);
      ctrl_logy_ss[j] = steps > 1 ? dot_self(ly - ctrl_logy_mean[j]) : 0.0;
      ctrl_ll_const[j] =
        -0.5 * (steps - 1) * log(2 * pi()) - 0.5 * log(steps) - sum(ly);

      nddc += steps;
    }
  }
}

parameters {
  // --------------------------------------------------------------------------
  // Population-level parameters
  // --------------------------------------------------------------------------

  real log_A0_pop;
  real log_gB_pop;

  // Weight on short component, on logit scale.
  real rho_pop;

  real log_r;

  // Ordered secretion rates:
  //   log_c_pop[1] = log_c_long
  //   log_c_pop[2] = log_c_short
  ordered[2] log_c_pop;

  // --------------------------------------------------------------------------
  // Baseline and case random-effect SDs
  // --------------------------------------------------------------------------

  real<lower=0> A0_sigma;
  real<lower=0> gB_sigma;

  real<lower=0> c_long_sigma;

  real<lower=0> rho_sigma;

  real<lower=0> c_short_sigma;

  // --------------------------------------------------------------------------
  // Case random effects
  // --------------------------------------------------------------------------

  vector[np] u_A0;
  vector[np] u_gB;

  vector[np] u_c_long;

  vector[np] u_rho;

  vector[np] u_c_short;

  // --------------------------------------------------------------------------
  // Observation noise
  // --------------------------------------------------------------------------

  real log_sigma_ctrl;
  real<lower=0> sigma_case;

  // --------------------------------------------------------------------------
  // Estimated delays
  // --------------------------------------------------------------------------

  vector<lower=delta_lower, upper=delta_upper>[n_delta] delta_free;
}

transformed parameters {
  vector[N] ypred;

  real<lower=0> sigma_ctrl = exp(log_sigma_ctrl);

  {
    vector[np] A0 = exp(log_A0_pop + u_A0 * A0_sigma);
    vector[np] gB = exp(log_gB_pop + u_gB * gB_sigma);

    vector[np] c_long = exp(log_c_pop[1] + u_c_long * c_long_sigma);
    vector[np] c_short = exp(log_c_pop[2] + u_c_short * c_short_sigma);

    real r = exp(log_r);

    vector[np] rho = inv_logit(rho_pop + u_rho * rho_sigma);

    int ndd = 1;

    for (i in 1:np) {
      int steps = nd[i];
      real delta_i = delta_id[i] > 0 ? delta_free[delta_id[i]] : 0.0;

      if (steps > 0) {
        ypred[ndd:(ndd + steps - 1)] =
          AbKin(
            t[ndd:(ndd + steps - 1)],
            A0[i],
            gB[i],
            delta_i,
            rho[i],
            c_long[i],
            r,
            c_short[i]
          );
      }

      ndd += steps;
    }
  }
}

model {
  // --------------------------------------------------------------------------
  // Scale priors, data-adaptive for nMFI
  // --------------------------------------------------------------------------

  log_A0_pop ~ normal(A0_prior_loc, A0_prior_sd);
  log_gB_pop ~ normal(log_gB_prior_loc, log_gB_prior_sd);

  // --------------------------------------------------------------------------
  // Population-level weight
  // --------------------------------------------------------------------------

  rho_pop ~ normal(rho_prior_loc, rho_prior_sd);

  // --------------------------------------------------------------------------
  // Rates, months^-1
  // --------------------------------------------------------------------------

  log_r ~ normal(log_r_prior_loc, log_r_prior_sd);

  log_c_pop[1] ~ normal(log_c_long_prior_loc, log_c_long_prior_sd);
  log_c_pop[2] ~ normal(log_c_short_prior_loc, log_c_short_prior_sd);

  // --------------------------------------------------------------------------
  // Random-effect SDs (half-normal)
  // --------------------------------------------------------------------------

  A0_sigma ~ normal(0, A0_sigma_prior_sd);
  gB_sigma ~ normal(0, gB_sigma_prior_sd);

  c_long_sigma ~ normal(0, c_long_sigma_prior_sd);
  c_short_sigma ~ normal(0, c_short_sigma_prior_sd);

  rho_sigma ~ normal(0, rho_sigma_prior_sd);

  // --------------------------------------------------------------------------
  // Standard normal random effects
  // --------------------------------------------------------------------------

  u_A0 ~ std_normal();
  u_gB ~ std_normal();

  u_c_long ~ student_t(u_c_long_prior_df, 0, 1);

  u_rho ~ std_normal();

  u_c_short ~ std_normal();

  // --------------------------------------------------------------------------
  // Observation-noise priors
  // --------------------------------------------------------------------------

  log_sigma_ctrl ~ normal(sigma_ctrl_prior_loc, sigma_ctrl_prior_sd);
  sigma_case ~ normal(0, sigma_case_prior_sd);

  // --------------------------------------------------------------------------
  // Delays
  // --------------------------------------------------------------------------

  delta_free ~ normal(0, 10 / 30.0);

  // --------------------------------------------------------------------------
  // Likelihood
  // --------------------------------------------------------------------------

  if (prior_only == 0) {
    if (n_obs > 0) {
      target += lognormal_lpdf(
        obs[idx_obs] | log(ypred[idx_obs]), sigma_case
      );
    }

    if (n_cens > 0) {
      target += lognormal_lccdf(
        upper_limit[idx_cens] | log(ypred[idx_cens]), sigma_case
      );
    }

    // Controls, baselines integrated out. For control j with n_j log-obs:
    //   mean(log y_j) ~ normal(log_A0_pop, sqrt(A0_sigma^2 + sigma_ctrl^2 / n_j))
    //   within-control SS contributes -(n_j - 1) log(sigma_ctrl) - SS / (2 sigma_ctrl^2)
    // Data-only constants are dropped here.
    target += normal_lpdf(
      ctrl_logy_mean |
      log_A0_pop,
      sqrt(square(A0_sigma) + square(sigma_ctrl) ./ ctrl_n)
    );
    target += -sum(ctrl_n - 1) * log_sigma_ctrl
              - sum(ctrl_logy_ss) / (2 * square(sigma_ctrl));
  }
}

generated quantities {
  vector[N] yrep_raw;
  vector[N] yrep;

  vector[nc] yrep_ctrl;

  vector[N] log_lik_case;
  // One log-lik value per control individual (identical to per-observation
  // when every control has a single sample).
  vector[npc] log_lik_ctrl;
  vector[N + npc] log_lik;

  vector[np] delta;

  vector[np] log_r_individual;
  vector[np] r_individual;

  vector[np] log_c_long_individual;
  vector[np] c_long_individual;

  vector[np] log_c_short_individual;
  vector[np] c_short_individual;

  vector[npc] u_A0_ctrl;
  vector[npc] log_A0_ctrl_individual;
  vector[npc] A0_ctrl_individual;

  vector[np] log_A0_individual;
  vector[np] A0_individual;
  vector[np] log_gB_individual;
  vector[np] gB_individual;

  vector[np] rho_individual;

  vector[np] w_short_individual;
  vector[np] w_long_individual;

  // Share of total antibody response (AUC) per component.
  vector[np] pi_short_individual;
  vector[np] pi_long_individual;
  real pi_short_pop;
  real pi_long_pop;

  for (i in 1:np) {
    delta[i] = delta_id[i] > 0 ? delta_free[delta_id[i]] : 0.0;

    log_r_individual[i] = log_r;
    r_individual[i] = exp(log_r);

    log_c_long_individual[i] = log_c_pop[1] + c_long_sigma * u_c_long[i];
    c_long_individual[i] = exp(log_c_long_individual[i]);

    log_c_short_individual[i] = log_c_pop[2] + c_short_sigma * u_c_short[i];
    c_short_individual[i] = exp(log_c_short_individual[i]);

    rho_individual[i] = inv_logit(rho_pop + rho_sigma * u_rho[i]);

    w_short_individual[i] = rho_individual[i];
    w_long_individual[i] = 1 - rho_individual[i];

    {
      real a_s = w_short_individual[i] / c_short_individual[i];
      real a_l = w_long_individual[i] / c_long_individual[i];

      pi_short_individual[i] = a_s / (a_s + a_l);
      pi_long_individual[i] = a_l / (a_s + a_l);
    }

    log_A0_individual[i] = log_A0_pop + A0_sigma * u_A0[i];
    A0_individual[i] = exp(log_A0_individual[i]);

    log_gB_individual[i] = log_gB_pop + gB_sigma * u_gB[i];
    gB_individual[i] = exp(log_gB_individual[i]);
  }

  {
    real rho_p = inv_logit(rho_pop);
    real a_s = rho_p / exp(log_c_pop[2]);
    real a_l = (1 - rho_p) / exp(log_c_pop[1]);

    pi_short_pop = a_s / (a_s + a_l);
    pi_long_pop = a_l / (a_s + a_l);
  }

  // Control baselines: draw from their conditional posterior given the
  // population parameters (normal-normal conjugacy on the log scale).
  {
    real tau2 = square(A0_sigma);
    real s2 = square(sigma_ctrl);
    int nddc = 1;

    for (j in 1:npc) {
      int steps = ndc[j];
      // With prior_only = 1 the controls carry no data, so draw from the prior.
      real prec_data = prior_only == 0 ? ctrl_n[j] / s2 : 0.0;
      real v = 1 / (1 / tau2 + prec_data);
      real m = v * (log_A0_pop / tau2 + prec_data * ctrl_logy_mean[j]);

      log_A0_ctrl_individual[j] = normal_rng(m, sqrt(v));
      A0_ctrl_individual[j] = exp(log_A0_ctrl_individual[j]);
      u_A0_ctrl[j] = (log_A0_ctrl_individual[j] - log_A0_pop) / A0_sigma;

      for (n in nddc:(nddc + steps - 1)) {
        yrep_ctrl[n] = lognormal_rng(log_A0_ctrl_individual[j], sigma_ctrl);
      }

      // Marginal log-likelihood of control j (all constants included).
      log_lik_ctrl[j] =
        ctrl_ll_const[j]
        + normal_lpdf(ctrl_logy_mean[j] |
                      log_A0_pop, sqrt(tau2 + s2 / ctrl_n[j]))
        - (ctrl_n[j] - 1) * log_sigma_ctrl
        - ctrl_logy_ss[j] / (2 * s2);

      log_lik[N + j] = log_lik_ctrl[j];

      nddc += steps;
    }
  }

  for (n in 1:N) {
    yrep_raw[n] = lognormal_rng(log(ypred[n]), sigma_case);

    if (is_upper_censored[n] == 1) {
      yrep[n] = fmin(yrep_raw[n], upper_limit[n]);
      log_lik_case[n] =
        lognormal_lccdf(upper_limit[n] | log(ypred[n]), sigma_case);
    } else {
      yrep[n] = yrep_raw[n];
      log_lik_case[n] =
        lognormal_lpdf(obs[n] | log(ypred[n]), sigma_case);
    }

    log_lik[n] = log_lik_case[n];
  }
}
