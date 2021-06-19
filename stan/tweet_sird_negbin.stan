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
 
      print("t=", t);
      print("SIRD=", y);
      print("beta, gamma, deathRate", theta);
      print("dS_dt, dI_Dt, dR_dt, dD_dt=", {dS_dt, dI_dt, dR_dt, dD_dt},
            " sum=", sum({dS_dt, dI_dt, dR_dt, dD_dt}));
  
      return {dS_dt, dI_dt, dR_dt, dD_dt};
  }

  real[,] my_sird(int days, real[] y, real[] theta, 
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
    /*
    print("beta, gamma, deathRate=", theta);  
    print("0 SIRD=", y," sum=", sum(y));
    */
    for (i in 1:days) {
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I ;
      real dR_dt =  gamma * I - deathRate * R;
      real dD_dt =  deathRate * R; 
    /*
    print("dS_dt, dI_Dt, dR_dt, dD_dt=", {dS_dt, dI_dt, dR_dt, dD_dt},
            " sum=", sum({dS_dt, dI_dt, dR_dt, dD_dt}));
            */
      S = dS_dt + S;
      I = dI_dt + I;
      R = dR_dt + R;
      D = dD_dt + D;
      day_counts[i] = {S, I, R, D};
/*
      print(i," SIRD=", day_counts[i], " sum=", sum(day_counts[i]));
      */
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
  int<lower = 0, upper = 1> run_SIR;
}
transformed data {
  real x_r[0]; //need for ODE function
  int x_i[1] = { N }; //need for ODE function

  //compartmentDays
  //sd
  //mean
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
  real<lower=.001> phi_inv;
  real<lower=.001> lambda_twitter;
  real<lower=.001> phi_twitter_inv;
}
transformed parameters{
  real y[n_days, n_compartments];
  real phi = 1. / phi_inv;
  real phi_twitter = 1. / phi_twitter_inv;
  matrix[n_days, n_compartments] daily_counts_ODE;
  real theta[3];
  theta[1] = beta;
  theta[2] = gamma;
  theta[3] = deaths;
//  y = integrate_ode_rk45(sird, y0, t0, ts, theta, x_r, x_i);
  y = my_sird(n_days, y0, theta, x_r, x_i);
  daily_counts_ODE = to_matrix(y);
  //if (check_ODE ==1 && sum(daily_counts_ODE) - n_days * N < 2) {
  //   
  
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
  sigma_compartment_noise ~ normal(0,20);
  sigma_i_compartment_noise ~ normal(0,20);
  sigma_twitter_noise ~ normal(0,20);
  lambda_twitter ~ normal(0,1);
  phi_twitter_inv ~ exponential(5);
  if (compute_likelihood == 1){ 
    for (i in 1:n_days) {
      if (run_SIR == 1) {
//        compartmentDays[i,compartment] ~ neg_binomial_2(daily_counts_ODE[i, compartment], phi);
        compartmentDays[i, compartment] ~ normal(daily_counts_ODE[i, compartment], 
                                       sigma_i_compartment_noise);
        /*compartmentDays[i,compartment] ~ normal(daily_counts_ODE[i, compartment], 
                                                sigma_compartment_noise);
                                                */
      }
      if (run_twitter == 1) {
        compartmentDays[i,tweetIndex] ~ normal(lambda_twitter * 
                                          daily_counts_ODE[i, tweetSourceIndex],
                                          sigma_twitter_noise);
        
        /*neg_binomial_2(lambda_twitter *
                                              daily_counts_ODE[i,tweetSourceIndex], 
                                              phi_twitter);*/
      }
    }
  }
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  real pred_tweets[n_days];
  real pred_i[n_days];
  matrix[n_days, n_compartments] ode_states = daily_counts_ODE;
  for (i in 1:n_days) {
      if (run_twitter == 1) {
          pred_tweets[i] = normal_rng(lambda_twitter *
                                   daily_counts_ODE[i, tweetSourceIndex],
                                   sigma_twitter_noise);
          
          /*
          neg_binomial_2_rng(lambda_twitter *
                                   daily_counts_ODE[i, tweetSourceIndex], 
                                   phi_twitter);
                                   */
      }
      if (run_SIR == 1) {
          pred_cases[i] = normal_rng(daily_counts_ODE[i, compartment], 
                                     sigma_compartment_noise);
                // neg_binomial_2_rng(daily_counts_ODE[i, compartment], 
                //                            phi);
         /* pred_i[i] = normal_rng(daily_counts_ODE[i,2], 
                      sigma_i_compartment_noise);
                      */
      }
  }
}
