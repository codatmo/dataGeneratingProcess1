/*
Modified from https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

*/

functions {
  real[] sird(real t, real[] y, real[] theta,
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real D = y[4];
      real N = x_i[1];

      real beta = theta[1];
      real gamma = theta[2];
      real deathRate = theta[3];

      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      real dD_dt =  deathRate * R;
  /*
            print("t=", t);

      print("SIRD=", y);
      print("N=", N);
      print("dS_dt, dI_Dt, dR_dt, dD_dt=", {dS_dt, dI_dt, dR_dt, dD_dt});
    */
      return {dS_dt, dI_dt, dR_dt, dD_dt};
  }
}
data {
  int<lower=1> n_days;
  real y0[4];
  real t0;
  real ts[n_days];
  int N;
  int death_count[n_days];
  int symptomaticTweets[n_days];
  int<lower = 0, upper = 1> compute_likelihood;
  int<lower = 0, upper = 1> use_twitter;
}
transformed data {
  real x_r[0]; //need for ODE function
  int x_i[1] = { N }; //need for ODE function
  int s_index = 1;
  int i_index = 2;
  int r_index = 3;
  int d_index = 4;
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> deaths;
  real<lower=.001> phi_inv;
  real<lower=.001> phi_twitter_inv;
  //real<lower=0> lag_twitter;
}
transformed parameters{
  real phi = 1. / phi_inv;
  real phi_twitter = 1. / phi_twitter_inv;

  // States to be recovered
  real state_estimate[n_days, 4];
  vector[n_days] state_S;
  vector[n_days] state_I;
  vector[n_days] state_R;
  vector[n_days] state_D;
  // int state_I_int[n_days];

  // ODE Stuff
  real theta[3];
  theta[1] = beta;
  theta[2] = gamma;
  theta[3] = deaths;

  state_estimate = integrate_ode_rk45(sird, y0, t0, ts, theta, x_r, x_i);

  // Populate States
  state_S = to_vector(state_estimate[, 1]);
  state_I = to_vector(state_estimate[, 2]);
  state_R = to_vector(state_estimate[, 3]);
  state_D = to_vector(state_estimate[, 4]);
  // state_I_int = to_int(state_estimate(state_I));

/*    print("y0=", y0, "t0=", t0, "ts=", ts, "theta=", theta, "x_r=", x_r,
          "x_i=", x_i, "state_estimate=", state_estimate);
    y0=[9999,1,0]
    t0=0ts=[1,2,3,4,5,6,7]
    theta=[0.0128292,0.049914]
    x_r=[]
    x_i=[10000]
    state_estimate=[[9998.99,0.963593,0.0489998],[9998.98,0.928512,0.0962157],
    [9998.96,0.894708,0.141713],[9998.95,0.862134,0.185553],[9998.94,0.830746,0.227797],
    [9998.93,0.800501,0.268504],[9998.92,0.771358,0.307728]]
*/

}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  deaths ~ normal(0.1, 0.1);
  phi_inv ~ exponential(5);
  //lag_twitter ~ normal(22, 5);
  phi_twitter_inv ~ exponential(5);
  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        // state_I[i] ~ neg_binomial_2(symptomaticTweets[i], phi_twitter);
      } else{
        death_count[i] ~ neg_binomial_2(state_D[i], phi);
      }
    }
  }
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  // real pred_cases[n_days];
  // real pred_tweets[n_days];
  // matrix[n_days, n_compartments] ode_states = daily_counts_ODE;
  // for (i in 1:n_days) {
  //     if (run_twitter == 1) {
  //         pred_tweets[i] = neg_binomial_2_rng(lambda_twitter *
  //                                  daily_counts_ODE[i, tweetSourceIndex],
  //                                  phi_twitter);
  //     }
  //     if (run_SIR == 1) {
  //         pred_cases[i] = neg_binomial_2_rng(daily_counts_ODE[i, compartment],
  //                                            phi);
  //     }
  // }
}
