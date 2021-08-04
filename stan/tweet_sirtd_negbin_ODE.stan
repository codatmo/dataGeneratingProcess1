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
  real prior_beta_mean;
  real prior_beta_std;
  real prior_omega_mean;
  real prior_omega_std;
  real prior_dI_mean;
  real prior_dI_std;
  real prior_dT_mean;
  real prior_dT_std;
  real prior_twitter_lambda;

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
  real<lower=0> twitter_rate;
}
transformed parameters{
  real phi = 1.0 / phi_inv;
  real phi_twitter = 1.0 / phi_twitter_inv;

  // States to be recovered
  vector[5] daily_counts_ODE[n_days];

  daily_counts_ODE = ode_rk45_tol(sird, y0, t0, ts,
                                1e-6, 1e-6, 1000, // if you want custom tolerances
                                beta, omega, dI, dT);
}
model {
  //priors
  beta ~ normal(prior_beta_mean, prior_beta_std);
  omega ~ normal(prior_omega_mean, prior_omega_std);
  dI ~ normal(prior_dI_mean, prior_dI_std);
  dT ~ normal(prior_dT_mean, prior_dT_std);
  phi_inv ~ exponential(5);
  phi_twitter_inv ~ exponential(5);
  twitter_rate ~ exponential(prior_twitter_lambda);

  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(twitter_rate * daily_counts_ODE[i, 2],
                                              phi_twitter);
      }
      death_count[i] ~ neg_binomial_2(daily_counts_ODE[i, 5], phi);
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
  state_S = to_vector(daily_counts_ODE[, 1]);
  state_I = to_vector(daily_counts_ODE[, 2]);
  state_R = to_vector(daily_counts_ODE[, 3]);
  state_T = to_vector(daily_counts_ODE[, 4]);
  state_D = to_vector(daily_counts_ODE[, 5]);

  // Tweets
  vector[n_days] state_tweets;
  for (i in 1:n_days) {
    state_tweets[i] = twitter_rate * state_I[i];
  }

  real R0 = beta * dI;
  int pred_deaths[n_days];
  int pred_tweets[n_days];
  real gamma = dI; // returning dI as gamma for breck's
  for (i in 1:n_days) {
     if (compute_likelihood == 1) {
          pred_deaths[i] = neg_binomial_2_rng(state_D[i], phi);
      }
      if (use_twitter == 1) {
          pred_tweets[i] = neg_binomial_2_rng(twitter_rate *
                                   state_I[i], phi_twitter);
      }
      else {
        pred_tweets[i] = 0;
      }
  }
}
