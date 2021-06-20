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
      real dI_dt =  beta * I * S / N - gamma * I ;
      real dR_dt =  gamma * I - deathRate * R;
      real dD_dt =  deathRate * R; 
 
      return {dS_dt, dI_dt, dR_dt, dD_dt};
  }

  real[,] block_sird(int days, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
    real day_counts[days,4];
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real D = y[4];
    real N = x_i[1];
    real beta = theta[1];
    real gamma = theta[2];
    real deathRate = theta[3];
    for (i in 1:days) {
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I ;
      real dR_dt =  gamma * I - deathRate * R;
      real dD_dt =  deathRate * R; 
      S = dS_dt + S;
      I = dI_dt + I;
      R = dR_dt + R;
      D = dD_dt + D;
      day_counts[i] = {S, I, R, D};
    }
    return day_counts;
  }
}

data {
  int<lower=1> n_days;
  int nDataCols;
  int<lower=1> n_compartments;
  real y0[4];
  real t0;
  real ts[n_days];
  int N;
  int compartmentDays[n_days, nDataCols];
  int compartment;
  int tweetIndex;
  int tweetSourceIndex;
  int<lower = 0, upper = 1> compute_likelihood;
  int<lower = 0, upper = 1> run_twitter;
  int<lower = 0, upper = 1> run_block_ODE;
  int<lower = 0, upper = 1> run_rk45_ODE;
}
transformed data {
  real x_r[0]; //need for ODE function
  int x_i[1] = { N }; //need for ODE function
  matrix[n_days, n_compartments] compartmentDaysM = to_matrix(compartmentDays);
  print("compiled");  
  print("compartment=", compartment);
}
parameters {
  real<lower = 0, upper = 1> gamma;
  real<lower=0, upper = 2> beta;
  real<lower=0, upper = 1> deaths;
  real<lower = 0> sigma_compartment_noise;
  real<lower = 0> sigma_i_compartment_noise;
  real<lower = 0> sigma_twitter_noise;
  real<lower=.001> lambda_twitter;
}
transformed parameters{
  real y[n_days, n_compartments];
  matrix[n_days, n_compartments] daily_counts_ODE;
  real theta[3];
  theta[1] = beta;
  theta[2] = gamma;
  theta[3] = deaths;
  if (run_rk45_ODE == 1 && run_block_ODE == 1) {
    reject("cannot run both rk45 and block ODEs");
  }
  if (run_rk45_ODE ==1 ) {
    y = integrate_ode_rk45(sird, y0, t0, ts, theta, x_r, x_i);
  }
  if (run_block_ODE == 1) {
    y = block_sird(n_days, y0, theta, x_r, x_i);
  }
  daily_counts_ODE = to_matrix(y);
}

model {
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  deaths ~ normal(0.1, 0.1);
  sigma_compartment_noise ~ normal(0,20);
  sigma_i_compartment_noise ~ normal(0,20);
  sigma_twitter_noise ~ normal(0,20);
  lambda_twitter ~ normal(0,1);
  if (compute_likelihood == 1) { 
    for (i in 1:n_days) {
      compartmentDays[i, compartment] ~ normal(daily_counts_ODE[i, compartment], 
                                       sigma_i_compartment_noise);

      if (run_twitter == 1) {
        compartmentDays[i,tweetIndex] ~ normal(lambda_twitter * 
                                          daily_counts_ODE[i, tweetSourceIndex],
                                          sigma_twitter_noise);
      }
    }
  }
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  vector[n_days] actual_cases = col(compartmentDaysM, compartment);
  vector[n_days] actual_tweets = col(compartmentDaysM, tweetIndex);
  real pred_deaths[n_days];
  real pred_tweets[n_days];
  real pred_i[n_days];
  matrix[n_days, n_compartments] ode_states = daily_counts_ODE;
  for (i in 1:n_days) {
      if (run_twitter == 1) {
          pred_tweets[i] = normal_rng(lambda_twitter *
                                   daily_counts_ODE[i, tweetSourceIndex],
                                   sigma_twitter_noise);
      }
      else {
        pred_tweets[i] = -100000;
      }
      pred_deaths[i] = normal_rng(daily_counts_ODE[i, compartment], 
                                 sigma_compartment_noise);
  }
}

