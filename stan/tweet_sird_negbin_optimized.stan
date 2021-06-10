/*
Modified from https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

*/

functions { // ODE new interface see here: https://mc-stan.org/users/documentation/case-studies/convert_odes.html
  vector sird(real time,
              vector y,
              real beta,
              real gamma,
              real deathRate,
              real dI,
              real N) {

      // Unpack state
      real S = y[1];
      real I = y[2];
      real R = y[3];
      real D = y[4];

      // derived parameters
      real infection = beta * I * S / N;
      real recovery = gamma * I;
      real death = deathRate * I * 1/dI;

      // ODE System
      vector[4] dydt;
      dydt[1] = -infection;                   // dS
      dydt[2] = infection - recovery - death; // dI
      dydt[3] = recovery;                     // dR
      dydt[4] = death;                        // dT

      return dydt;
  }
}
data {
  int<lower=1> n_days;
  vector[4] y0;
  real t0;
  real ts[n_days];
  real N;
  int death_count[n_days];
  int symptomaticTweets[n_days];
  int<lower=0, upper=1> compute_likelihood;
  int<lower=0, upper=1> use_twitter;
}
transformed data {
  int s_index = 1;
  int i_index = 2;
  int r_index = 3;
  int d_index = 4;
}
parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> death_rate;
  real<lower=0> dI;
  real<lower=0.001> phi_inv;
  real<lower=0.001> phi_twitter_inv;
  real<lower=0, upper=1> proportion_twitter;
  //real<lower=0> lag_twitter;
  real<lower=0> I_noise;
  real<lower=0> twitter_noise;
}
transformed parameters{
  // I'm not giving up negbin yet
  real phi = 1.0 / phi_inv;
  real phi_twitter = 1.0 / phi_twitter_inv;

  // States to be recovered
  vector<lower=0>[4] state_estimate[n_days];
  vector[n_days] state_S;
  vector[n_days] state_I;
  vector[n_days] state_R;
  vector[n_days] state_D;
  // int state_I_int[n_days];

  // ODE RK45 (Stan's Implementation)
  // state_estimate = ode_rk45(sird, y0, t0, ts,
  //                           beta, gamma, death_rate, dI, N);
  state_estimate = ode_rk45_tol(sird, y0, t0, ts,
                                1e-6, 1e-6, 1000, // if you want custom tolerances
                                beta, gamma, death_rate, dI, N);

  // Populate States
  state_S = to_vector(state_estimate[, 1]);
  state_I = to_vector(state_estimate[, 2]);
  state_R = to_vector(state_estimate[, 3]);
  state_D = to_vector(state_estimate[, 4]);

  // We might need this in the future
  vector[n_days-1] daily_infections;
  vector[n_days-1] daily_deaths;
  daily_infections = state_S[:n_days-1] - state_S[2:] + 1E-4;
  daily_deaths = state_D[2:] - state_D[:n_days-1] + 1E-4; // I had problems with negativity
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  death_rate ~ normal(0.1, 0.1);
  dI ~ normal(7, 1); // gammaResRateSim = 1 / 7
  phi_inv ~ exponential(5);
  phi_twitter_inv ~ exponential(5);
  proportion_twitter ~ beta(1, 1); // Beta is a better prior for proportions
  //lag_twitter ~ normal(22, 5); // in the future I want to put a lag

  // Compartment Noises
  I_noise ~ normal(0, 20);       // 20 is too high...
  twitter_noise ~ normal(0, 20); // 20 is too high...

  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(state_I[i], phi_twitter); // and a incorporate lag
        // symptomaticTweets[i] ~ normal(proportion_twitter * state_I[i], twitter_noise);
      } else {
        death_count[i] ~ neg_binomial_2(state_D[i], phi);
        // death_count[i] ~ normal(state_D[i], I_noise);
      }
    }
  }
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  real pred_tweets[n_days];
  vector[4] ode_states[n_days]  = state_estimate;
  for (i in 1:n_days) {
     if (compute_likelihood == 1) {
          pred_cases[i] = normal_rng(state_D[i], twitter_noise);
      }
      if (use_twitter == 1) {
          pred_tweets[i] = normal_rng(proportion_twitter *
                                   state_I[i], I_noise);
      }
  }
}
