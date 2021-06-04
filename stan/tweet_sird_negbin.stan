/*
Modified from https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

*/

functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real D = y[4];
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      real deathRate = theta[3];
      // real alpha_tweet_generation = theta[3];
      
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      real dD_dt =  deathRate * R; 
      // real dT_dt = alpha_tweet_generation * ;I
      
      return {dS_dt, dI_dt, dR_dt, dD_dt};
  }
}
data {
  int<lower=1> n_days;
   real y0[4];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
  int symptomaticTweets[n_days];
  int<lower = 0, upper = 1> compute_likelihood;
  int<lower = 0, upper = 1> run_twitter;
  int<lower = 0, upper = 1> run_SIR;
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
  real<lower=.001> lambda_twitter;
  real<lower=.001> phi_twitter_inv;
}
transformed parameters{
  real y[n_days, 4];
  real phi = 1. / phi_inv;
  real phi_twitter = 1. / phi_twitter_inv;
  vector[n_days] infected_daily_counts_ODE;
    vector[n_days] deaths_daily_counts_ODE;
    real theta[3];
    theta[1] = beta;
    theta[2] = gamma;
    theta[3] = deaths;
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
    infected_daily_counts_ODE = col(to_matrix(y), i_index);
    deaths_daily_counts_ODE = col(to_matrix(y), d_index);
/*    print("y0=", y0, "t0=", t0, "ts=", ts, "theta=", theta, "x_r=", x_r, 
          "x_i=", x_i, "y=", y);
    y0=[9999,1,0]
    t0=0ts=[1,2,3,4,5,6,7]
    theta=[0.0128292,0.049914]
    x_r=[]
    x_i=[10000]
    y=[[9998.99,0.963593,0.0489998],[9998.98,0.928512,0.0962157],
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
  lambda_twitter ~ normal(0,1);
  phi_twitter_inv ~ exponential(5);
  if (compute_likelihood == 1){ 
    for (i in 1:n_days) {
      if (run_SIR == 1) {
        cases[i] ~ neg_binomial_2(deaths_daily_counts_ODE[i], phi);
      }
      if (run_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(lambda_twitter*
                                              infected_daily_counts_ODE[i], 
                                              phi_twitter);
      }
    }
  }
}
generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  real pred_tweets[n_days];
  for (i in 1:n_days) {
      if (run_twitter == 1) {
          pred_tweets[i] = neg_binomial_2_rng(lambda_twitter*infected_daily_counts_ODE[i], phi_twitter);
      }
      if (run_SIR == 1) {
          pred_cases[i] = neg_binomial_2_rng(infected_daily_counts_ODE[i], phi);
      }
  }
}
