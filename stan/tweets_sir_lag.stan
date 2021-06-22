functions {
  //version which uses bisectioning search
  int find_interval_elem(real x, vector sorted, int start_ind) {
    int res;
    int N = num_elements(sorted);
    int max_iter = 100 * N;
    int left_ind = start_ind;
    int right_ind = N;
    real left = sorted[left_ind ] - x;
    real right = sorted[right_ind] - x;
    int iter = 1;

    if(N == 0) return(0);
    if(0 == left) return(left_ind);
    if(0 < left) return(left_ind-1);
    if(0 >= right) return(right_ind);

    while((right_ind - left_ind) > 1  && iter != max_iter) {
      int mid_ind;
      real mid;
      // is there a controlled way without being yelled at with a
      // warning?
      mid_ind = (left_ind + right_ind) / 2;
      mid = sorted[mid_ind] - x;
      if (mid == 0) return(mid_ind-1);
      if (left  * mid < 0) { right = mid; right_ind = mid_ind; }
      if (right * mid < 0) { left  = mid; left_ind  = mid_ind; }
      iter = iter + 1;
    }
    if(iter == max_iter)
      print("Maximum number of iterations reached.");
    return(left_ind);
  }

  // ODE new interface see here: https://mc-stan.org/users/documentation/case-studies/convert_odes.html
  vector sird(real time,
              vector y,
              real beta,
              real gamma,
              real deathRate,
              real N) {

      // Unpack state
      real S = y[1];
      real I = y[2];
      real R = y[3];
      real D = y[4];

      // derived parameters
      real infection = beta * I * S / N;
      real recovery = gamma * I;
      real death = deathRate * I * gamma;

      // ODE System
      vector[4] dydt;
      dydt[1] = -infection;                   // dS
      dydt[2] = infection - recovery - death; // dI
      dydt[3] = recovery;                     // dR
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
  int max_lag = 7; // in real workd would be approx 23 + some wiggle room
}
parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> death_rate;
  real<lower=0.001> phi_inv;
  real<lower=0.001> phi_twitter_inv;
  real<lower=0, upper=1> proportion_twitter;
  simplex[max_lag+1] lag_weights_twitter_symptoms;
  // real<lower=0> I_noise;
  // real<lower=0> twitter_noise;
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

  if (trapezoidal_solver) {
    state_estimate = ode_explicit_trapezoidal(y0, t0, ts,
                                              beta, gamma, death_rate, N);
  } else {
    // ODE RK45 (Stan's Implementation)
    // state_estimate = ode_rk45(sird, y0, t0, ts,
    //                           beta, gamma, death_rate, N);
    state_estimate = ode_rk45_tol(sird, y0, t0, ts,
                                  1e-6, 1e-6, 1000, // if you want custom tolerances
                                  beta, gamma, death_rate, N);
  }

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

  // Twitter lag scratch
  twitter_symptoms_lagged_daily_infections = lag_weights_twitter_symptoms[1] * daily_infections;
  for (i in 1:max_lag) {
    twitter_symptoms_lagged_daily_infections += lag_weights_twitter_symptoms[i+1]*
                                          append_row(rep_vector(0.0, i), daily_infections[:n_days-i]);
  }

  daily_twitter_symptoms = rep_vector(0.0, n_days);
  daily_twitter_symptoms = twitter_symptoms_lagged_daily_infections * proportion_twitter;
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  death_rate ~ normal(0.1, 0.1);
  phi_inv ~ exponential(5);
  phi_twitter_inv ~ exponential(5);
  proportion_twitter ~ beta(1, 1); // Beta is a better prior for proportions
  lag_weights_twitter_symptoms ~ dirichlet(rep_vector(0.1, max_lag+1));

  // Compartment Noises
  // I_noise ~ normal(0, 20);       // 20 is too high...
  // twitter_noise ~ normal(0, 20); // 20 is too high...

  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(daily_twitter_symptoms[i], phi_twitter);
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
  vector[4] ode_states[n_days]  = state_estimate;
  for (i in 1:n_days) {
     if (compute_likelihood == 1) {
          // pred_deaths[i] = normal_rng(state_D[i], I_noise);
          pred_deaths[i] = neg_binomial_2_rng(state_D[i], phi);
      }
      if (use_twitter == 1) {
          // pred_tweets[i] = normal_rng(proportion_twitter *
          //                          state_I[i], I_noise);
          pred_tweets[i] = neg_binomial_2_rng(daily_twitter_symptoms[i],
                                              phi_twitter);
      }
  }
}
