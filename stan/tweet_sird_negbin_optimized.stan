functions { // ODE new interface see here: https://mc-stan.org/users/documentation/case-studies/convert_odes.html
  vector sird(real time,
              vector y,
              real beta,
              real omega,
              real dI,
              real dT) {

      // Unpack state
      real S = y[1];
      real I = y[2];
      real R = y[3];
      real T = y[4];
      real D = y[5];
      real N = S+I+R+T+D;

      // ODE System
      vector[5] dydt;
      dydt[1] = -beta * I * S / N;            // dS
      dydt[2] = beta * I * S / N - 1/dI * I;  // dI
      dydt[3] = 1/dI * I * (1 - omega);       // dR
      dydt[4] = 1/dI * I * omega - 1/dT * T;  // dT
      dydt[5] = 1/dT * T;                     // dD

      return dydt;
  }
}
data {
  int<lower=1> n_days;
  vector[5] y0;
  real t0;
  real ts[n_days];
  int death_count[n_days];
  int symptomaticTweets[n_days];
  int<lower=0, upper=1> compute_likelihood;
  int<lower=0, upper=1> use_twitter;
}
transformed data {
  // if necessary
}
parameters {
  // ODE Stuff
  real<lower=0> beta;
  real<lower=0> omega;
  real<lower=0> dI;
  real<lower=0> dT;

  real<lower=0.001> phi_inv;
  real<lower=0.001> phi_twitter_inv;
  real<lower=0, upper=1> proportion_twitter;
}
transformed parameters{
  real phi = 1.0 / phi_inv;
  real phi_twitter = 1.0 / phi_twitter_inv;

  // States to be recovered
  vector[5] state_estimate[n_days];

  state_estimate = ode_rk45_tol(sird, y0, t0, ts,
                                1e-6, 1e-6, 1000, // if you want custom tolerances
                                beta, omega, dI, dT);
}
model {
  //priors
  beta ~ normal(2, 1);
  omega ~ normal(0.4, 0.5);
  dI ~ normal(7, 2);
  dT ~ normal(10, 2);
  phi_inv ~ exponential(5);
  phi_twitter_inv ~ exponential(5);
  proportion_twitter ~ beta(1, 2); // Beta is a better prior for proportions

  // Compartment Noises
  I_noise ~ normal(0, 20);       // 20 is too high...
  twitter_noise ~ normal(0, 20); // 20 is too high...

  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(proportion_twitter * state_estimate[i, 2],
                                              phi_twitter);
        // symptomaticTweets ~ normal(proportion_twitter * state_estimate[i, 2], twitter_noise);
      }
      death_count[i] ~ neg_binomial_2(state_estimate[i, 5], phi);
      // death_count[i] ~ normal(state_estimate[i, 5], I_noise);
    }
  }
}

generated quantities {
  // States to be recovered
  vector[n_days] state_S;
  vector[n_days] state_I;
  vector[n_days] state_R;
  vector[n_days] state_T;
  vector[n_days] state_D;

  // Populate States
  state_S = to_vector(state_estimate[, 1]);
  state_I = to_vector(state_estimate[, 2]);
  state_R = to_vector(state_estimate[, 3]);
  state_T = to_vector(state_estimate[, 4]);
  state_D = to_vector(state_estimate[, 5]);

  // Tweets
  vector[n_days] state_tweets;
  for (i in 1:n_days) {
    state_tweets[i] = proportion_twitter * state_I[i];
  }

  real R0 = beta * dI;
  int pred_deaths[n_days];
  int pred_tweets[n_days];
  for (i in 1:n_days) {
     if (compute_likelihood == 1) {
          pred_deaths[i] = neg_binomial_2_rng(state_D[i], phi);
      }
      if (use_twitter == 1) {
          pred_tweets[i] = neg_binomial_2_rng(proportion_twitter *
                                   state_I[i], phi_twitter);
      }
  }
}
