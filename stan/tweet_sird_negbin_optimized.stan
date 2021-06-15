functions { // ODE new interface see here: https://mc-stan.org/users/documentation/case-studies/convert_odes.html
  vector sird(real time,
              vector y,
              real beta,
              real gamma,
              real omega,
              real N) {

      // Unpack state
      real S = y[1];
      real I = y[2];
      real R = y[3];
      real D = y[4];

      // derived parameters
      real infection = beta * I * S / N;
      real recovery = gamma * I;
      real death = omega * R;

      // ODE System
      vector[4] dydt;
      dydt[1] = -infection;                   // dS
      dydt[2] = infection - recovery;         // dI
      dydt[3] = recovery - death;             // dR
      dydt[4] = death;                        // dT

      return dydt;
  }

  vector[] ode_explicit_trapezoidal(vector initial_state,
                                            real initial_time,
                                            real[] times,
                                            real beta,
                                            real gamma,
                                            real deathRate,
                                            real N) {
    real h;
    vector[size(initial_state)] dstate_dt_initial_time;
    vector[size(initial_state)] dstate_dt_tidx;
    vector[size(initial_state)] k;
    vector[size(initial_state)] state_estimate[size(times)];

    h = times[1] - initial_time;
    dstate_dt_initial_time = sird(initial_time, initial_state, beta, gamma, deathRate, N);
    k = h * dstate_dt_initial_time;
    state_estimate[1,] = initial_state + h * (dstate_dt_initial_time + sird(times[1], initial_state + k, beta, gamma, deathRate, N)) / 2;

    for (tidx in 1:size(times)-1) {
      h = (times[tidx+1] - times[tidx]);
      dstate_dt_tidx = sird(times[tidx], state_estimate[tidx], beta, gamma, deathRate, N);
      k = h * dstate_dt_tidx;
      state_estimate[tidx+1,] = state_estimate[tidx,] + h * (dstate_dt_tidx + sird(times[tidx+1], state_estimate[tidx,] + k, beta, gamma, deathRate, N))/2;
    }

    return state_estimate;
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
  int<lower=0, upper=1> trapezoidal_solver;
}
transformed data {

}
parameters {
  real<lower=0.001> beta;
  real<lower=0.001> gamma;
  real<lower=0.001> omega;
  real<lower=0.001> phi_inv;
  real<lower=0.001> phi_twitter_inv;
  real<lower=0, upper=1> proportion_twitter;

  real<lower=0> I_noise;
  real<lower=0> twitter_noise;
}
transformed parameters{
  real phi = 1.0 / phi_inv;
  real phi_twitter = 1.0 / phi_twitter_inv;

  // States to be recovered
  vector<lower=0>[4] state_estimate[n_days];
  vector[n_days] state_S;
  vector[n_days] state_I;
  vector[n_days] state_R;
  vector[n_days] state_D;

  if (trapezoidal_solver) {
    state_estimate = ode_explicit_trapezoidal(y0, t0, ts,
                                              beta, gamma, omega, N);
  } else {
    // ODE RK45 (Stan's Implementation)
    // state_estimate = ode_rk45(sird, y0, t0, ts,
    //                           beta, gamma, omega, N);
    state_estimate = ode_rk45_tol(sird, y0, t0, ts,
                                  1e-6, 1e-6, 1000, // if you want custom tolerances
                                  beta, gamma, omega, N);
  }

  // Populate States
  state_S = to_vector(state_estimate[, 1]);
  state_I = to_vector(state_estimate[, 2]);
  state_R = to_vector(state_estimate[, 3]);
  state_D = to_vector(state_estimate[, 4]);

  // Daily Stuff
  vector[n_days] daily_infections;
  vector[n_days] daily_deaths;
  daily_infections = append_row(y0[2], to_vector(state_S[:n_days-1] - state_S[2:] + 1E-4));
  daily_deaths = append_row(y0[4], to_vector(state_D[2:] - state_D[:n_days-1] + 1E-4)); // I had problems with negativity
}
model {
  //priors
  beta ~ normal(0.3, 0.1);
  gamma ~ normal(0.14, 0.1);
  omega ~ normal(0.1, 0.1);
  phi_inv ~ exponential(5);
  phi_twitter_inv ~ exponential(5);
  proportion_twitter ~ beta(1, 1); // Beta is a better prior for proportions

  // Compartment Noises
  I_noise ~ normal(0, 20);       // 20 is too high...
  twitter_noise ~ normal(0, 20); // 20 is too high...

  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(proportion_twitter * state_I[i],
                                              phi_twitter);
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
