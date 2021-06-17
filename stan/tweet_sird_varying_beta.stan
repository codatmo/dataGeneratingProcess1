functions { // ODE new interface see here: https://mc-stan.org/users/documentation/case-studies/convert_odes.html
  real[] sird(real time,
              real[] state,
              real[] params,
              real[] real_data,
              int[] integer_data) {

      // Unpack integer data values
    int T = integer_data[1];
    int n_beta_pieces = integer_data[2];

    // Unpack real data values
    real beta_left_t[n_beta_pieces] = real_data[1:n_beta_pieces];
    real beta_right_t[n_beta_pieces] = real_data[n_beta_pieces+1:2*n_beta_pieces];
    real population = real_data[2*n_beta_pieces+1];

    // Unpack parameter values
    real beta_left[n_beta_pieces] = params[1:n_beta_pieces];
    real grad_beta[n_beta_pieces] = params[n_beta_pieces+1:2*n_beta_pieces];
    real gamma = params[2*n_beta_pieces+1];
    real kappa = params[2*n_beta_pieces+2];
    real omega = params[2*n_beta_pieces+3];

    // Unpack state
    real S = state[1];
    real I = state[2];
    real T_ = state[3];
    real D = state[4];

    real infection_rate;
    real gammaI = gamma * I;
    real kappaT = kappa * T_;

    real dS_dt;
    real dI_dt;
    real dR_dt;
    real dT_dt;
    real dD_dt;

    for (i in 1:n_beta_pieces) {
      if(time >= beta_left_t[i] && time < beta_right_t[i]) {
        real beta = grad_beta[i] * (time - beta_left_t[i]) + beta_left[i];
        infection_rate = beta * I * S / population;
      }
    }

    dS_dt = -infection_rate;
    dI_dt = infection_rate - gammaI;
    dT_dt = gammaI * omega - kappaT;
    dD_dt = kappaT;

    return {dS_dt, dI_dt, dT_dt, dD_dt};
  }

  real[ , ] integrate_ode_explicit_trapezoidal(real[] initial_state, real initial_time, real[] times, real[] params, real[] real_data, int[] integer_data) {
    real h;
    vector[size(initial_state)] dstate_dt_initial_time;
    vector[size(initial_state)] dstate_dt_tidx;
    vector[size(initial_state)] k;
    real state_estimate[size(times),size(initial_state)];

    h = times[1] - initial_time;
    dstate_dt_initial_time = to_vector(sird(initial_time, initial_state, params, real_data, integer_data));
    k = h*dstate_dt_initial_time;
    state_estimate[1,] = to_array_1d(to_vector(initial_state) + h*(dstate_dt_initial_time + to_vector(sird(times[1], to_array_1d(to_vector(initial_state)+k), params, real_data, integer_data)))/2);

    for (tidx in 1:size(times)-1) {
      h = (times[tidx+1] - times[tidx]);
      dstate_dt_tidx = to_vector(sird(times[tidx], state_estimate[tidx], params, real_data, integer_data));
      k = h*dstate_dt_tidx;
      state_estimate[tidx+1,] = to_array_1d(to_vector(state_estimate[tidx,]) + h*(dstate_dt_tidx + to_vector(sird(times[tidx+1], to_array_1d(to_vector(state_estimate[tidx,])+k), params, real_data, integer_data)))/2);
    }

    return state_estimate;
  }
}
data {
  int<lower=1> n_days;
  vector[4] y0;
  real t0;
  real ts[n_days];
  int death_count[n_days];
  int symptomaticTweets[n_days];
  int<lower=1> n_beta_pieces;
  real<lower=0> beta_left_t[n_beta_pieces];
  real<lower=0> beta_right_t[n_beta_pieces];
  int real_data_length;
  real real_data[real_data_length];
  int integer_data_length;
  int integer_data[integer_data_length];
  int<lower=0, upper=1> compute_likelihood;
  int<lower=0, upper=1> use_twitter;
}
transformed data {
  // if necessary
}
parameters {
  // ODE Stuff
  real<lower=0> beta_left[n_beta_pieces];
  real<lower=0> beta_right[n_beta_pieces];
  real<lower=0> omega;
  real<lower=0> dI;
  real<lower=0> dT;

  real<lower=0.001> phi_inv;
  real<lower=0.001> phi_twitter_inv;
  real<lower=0, upper=1> proportion_twitter;

  real<lower=0> I_noise;
  real<lower=0> twitter_noise;
}
transformed parameters{
  real grad_beta[n_beta_pieces];
  real gamma;
  real kappa;
  real phi = 1.0 / phi_inv;
  real phi_twitter = 1.0 / phi_twitter_inv;

  real state_estimate[n_days, 4];
  vector[n_days] S;
  vector[n_days] I;
  vector[n_days] T;
  vector[n_days] D;

  grad_beta = to_array_1d((to_vector(beta_right) - to_vector(beta_left))./(to_vector(beta_right_t) -
              to_vector(beta_left_t)));
  gamma = 1.0/dI;
  kappa = 1.0/dT;

  // States to be recovered
  {
    real params[2*n_beta_pieces+4];
    params[1:n_beta_pieces] = beta_left;
    params[n_beta_pieces+1:2*n_beta_pieces] = grad_beta;
    params[2*n_beta_pieces+1] = gamma;
    params[2*n_beta_pieces+2] = kappa;
    params[2*n_beta_pieces+3] = omega;

    state_estimate = integrate_ode_explicit_trapezoidal(to_array_1d(y0), t0, ts, params, real_data, integer_data);
  }
  S = to_vector(state_estimate[, 1]);
  I = to_vector(state_estimate[, 2]);
  T = to_vector(state_estimate[, 3]);
  D = to_vector(state_estimate[, 4]);
}
model {
  //priors
  beta_left ~ normal(0, 0.5);
  beta_right ~ normal(0, 0.5);
  omega ~ beta(1, 10);
  dI ~ normal(7, 1);
  dT ~ normal(10, 1);
  phi_inv ~ exponential(5);
  phi_twitter_inv ~ exponential(5);
  proportion_twitter ~ beta(1, 1); // Beta is a better prior for proportions

  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(proportion_twitter * I[i],
                                              phi_twitter);
      }
      death_count[i] ~ neg_binomial_2(D[i], phi);
    }
  }
}

generated quantities {
  int pred_deaths[n_days];
  int pred_tweets[n_days];
  for (i in 1:n_days) {
     if (compute_likelihood == 1) {
          pred_deaths[i] = neg_binomial_2_rng(D[i], phi);
      }
      if (use_twitter == 1) {
          pred_tweets[i] = neg_binomial_2_rng(proportion_twitter * I[i], phi_twitter);
      }
  }
}
