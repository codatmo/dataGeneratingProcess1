/*
Modified from https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html
*/

// #include /home/breck/git/codatmo/dataGeneratingProcess1/stan/ode_solvers.stan

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
             real[] x_r, int[] x_i, int debug) {
    real day_counts[days,4];
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real D = y[4];
    real N = x_i[1];
    real beta = theta[1];
    real gamma = theta[2];
    real deathRate = theta[3];
    if (debug == 1) {
      print("beta, gamma, deathRate=", theta);  
      print("0 SIRD=", y," sum=", sum(y));
    }  
    for (i in 1:days) {
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I ;
      real dR_dt =  gamma * I - deathRate * R;
      real dD_dt =  deathRate * R; 
    if (debug ==1) {  
      print("dS_dt, dI_Dt, dR_dt, dD_dt=", {dS_dt, dI_dt, dR_dt, dD_dt},
            " sum=", sum({dS_dt, dI_dt, dR_dt, dD_dt}));
    }          
      S = dS_dt + S;
      I = dI_dt + I;
      R = dR_dt + R;
      D = dD_dt + D;
      day_counts[i] = {S, I, R, D};
     if (debug == 1) { 
        print(i," SIRD=", day_counts[i], " sum=", sum(day_counts[i]));
     }
    }
    return day_counts;
  }

}
data {
  int<lower=1> n_days;
  real sDay1;
  real iDay1;
  real rDay1;
  real dDay1;
  int NPop;
  vector[n_days] tweets;
  vector[n_days] deaths;
  int<lower = 0, upper = 1> compute_likelihood;
  int<lower = 0, upper = 1> run_twitter;
  int<lower = 0, upper = 1> run_block_ODE;
  int<lower = 0, upper = 1> run_rk45_ODE;
}
transformed data {
  real x_r[0]; //need for ODE function
  int x_i[1] = { NPop }; //need for ODE function
  int debug = 0;
  int debug2 = 0;
  int n_compartments = 4;
  int sCompartment = 1;
  int iCompartment = 2;
  int rCompartment = 3;
  int dCompartment = 4;
  real sdTweets = sd(tweets);
  real meanTweets = mean(tweets);
  vector[n_days] tweets_s_c = (tweets)/sdTweets;
  real sdDeaths = sd(deaths);
  real meanDeaths = mean(deaths);
  vector[n_days] deaths_s_c = (deaths)/sdDeaths;
  real ts[n_days];
  ts[1] = 1.0;
  for (i in 2:n_days) {
      ts[i] = ts[i - 1] + 1;
  }
  print("tweets_s_c", tweets_s_c);
  print("deaths_s_c", deaths_s_c);
}
parameters {
  real<lower = 0, upper = 1> gamma;
  real<lower=0, upper = 2> beta;
  real<lower=0, upper = 1> deathRate;
  real<lower = 0> sigma_sd_beta;
  real<lower = 0> sigma_sd_gamma;
  real<lower = 0> sigma_sd_deaths;
  real<lower = 0> sigma_compartment_noise;
  real<lower = 0> sigma_twitter_noise;
  real<lower=.001> lambda_twitter;
}
transformed parameters{
  real compartmentStartValues[4] = {sDay1, iDay1, rDay1, dDay1};  
  real y[n_days, 4];
  matrix[n_days, 4] daily_counts_ODE;
  real theta[3];
  theta[1] = beta;
  theta[2] = gamma;
  theta[3] = deathRate;
  
  if (run_rk45_ODE == 1 && run_block_ODE == 1) {
    reject("cannot run both rk45 and block ODEs");
  }
  if (run_rk45_ODE == 1 ) {
    y = integrate_ode_rk45(sird, compartmentStartValues , 0.0, ts, theta, x_r, x_i);
  }
  if (run_block_ODE == 1) {
    y = block_sird(n_days, compartmentStartValues, theta, x_r, x_i, debug);
  }
  daily_counts_ODE = to_matrix(y);
}

model {
  sigma_sd_beta ~ exponential(10);
  sigma_sd_gamma ~ exponential(10);
  sigma_sd_deaths ~ exponential(10);
  beta ~ normal(0, sigma_sd_beta);
  gamma ~ normal(0, sigma_sd_gamma);
  deathRate ~ normal(0, sigma_sd_deaths);
  lambda_twitter ~ normal(0,1);
  if (compute_likelihood == 1) { 
    for (i in 1:n_days) {
      	deaths_s_c[i] ~ normal(daily_counts_ODE[i, dCompartment], 
                                       sigma_compartment_noise);
      if (run_twitter == 1) {
        tweets_s_c[i] ~ normal(lambda_twitter * daily_counts_ODE[i, iCompartment],
                                          sigma_twitter_noise);
      }
    }
  }
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_deaths[n_days];
  real pred_tweets[n_days];
  real pred_i[n_days];
  matrix[n_days, n_compartments] ode_states = daily_counts_ODE;
  for (i in 1:n_days) {
      if (run_twitter == 1) {
          pred_tweets[i] = sdTweets * (normal_rng(lambda_twitter *
                                   daily_counts_ODE[i, dCompartment],
                                   sigma_twitter_noise));
      }
      else {
        pred_tweets[i] = 0;
      }
      pred_deaths[i] = sdDeaths * (normal_rng(daily_counts_ODE[i, iCompartment], 
                                 sigma_compartment_noise));
  }
}

