data {
  int<lower=1> n_days;
  vector[5] y0;
  real t0;
  real ts[n_days];
  int death_count[n_days];
  int symptomaticTweets[n_days];
  int<lower=0, upper=1> compute_likelihood;
  int<lower=0, upper=1> use_twitter;
  real<lower=0> beta_regularization;
  real prior_omega_mean;
  real prior_omega_std;
  real prior_dI_mean;
  real prior_dI_std;
  real prior_dT_mean;
  real prior_dT_std;
  real prior_twitter_lambda;
}
transformed data {
  real population = sum(y0);

  // Model needs new deaths
  int new_deaths[n_days];
  new_deaths[1] = 0;
  for (i in 2:n_days) {
    new_deaths[i] = death_count[i] - death_count[i-1];
  }
  // Varying Weekly Beta Stuff
  int n_weeks = n_days %/% 7 + min(1, n_days % 7);
  int new_weekly_deaths[n_weeks];
  for(week in 1:n_weeks){
    int start = 1+7*(week-1);
    int end = min(start+6, n_days);
    new_weekly_deaths[week] = sum(new_deaths[start:end]);
  }
}
parameters {
  /*
    The elements of a simplex sum to one.
    The first n_days entries represent daily new infections as a fraction
    of the total population while the last entry represents the proportion that
    remains susceptible at the end.
  */
  simplex[n_days+1] unit_dS;

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

  vector[n_days] daily_infections = population * unit_dS[:n_days];
  vector[n_days] daily_deaths;
  vector[n_weeks] weekly_deaths;
  vector[n_days] beta;
  vector[n_days] effective_reproduction_number;

  // States to be recovered
  vector[n_days] state_S;
  vector[n_days] state_I;
  vector[n_days] state_R;
  vector[n_days] state_T;
  vector[n_days] state_D;

  // matrix_exp
  if(compute_likelihood){
    vector[4] state = [
        y0[2],
        y0[3],
        y0[4],
        y0[5]
    ]';
    matrix[4, 4] transition_matrix = matrix_exp([
    //[I                  R, ,T     ,D]
      [-1/dI             ,0  ,0      ,0],//I
      [+1/dI*(1 - omega) ,0  ,0      ,0],//R
      [+1/dI*omega       ,0  ,-1/dT  ,0],//T
      [0                 ,0  ,+1/dT  ,0] //D
    ]);
    real S = population;
    real last_D;
    int weekday;
    for(i in 1:n_days){
      weekday = 1+(i-1) % 7;
      last_D = state[4];
      S -= daily_infections[i];
      state[1] += daily_infections[i];
      state = transition_matrix * state;
      daily_deaths[i] = state[4] - last_D;
      beta[i] = daily_infections[i] * population / (S * state[1]); // S * I
      effective_reproduction_number[i] = daily_infections[i] / state[1] * dI; // I

      // Populate States
      state_S[i] = S;
      state_I[i] = state[1];
      state_R[i] = state[2];
      state_T[i] = state[3];
      state_D[i] = state[4];
    }
    for(week in 1:n_weeks){
      int start = 1+7*(week-1);
      int end = min(start+6, n_days);
      weekly_deaths[week] = sum(daily_deaths[start:end]);
    }
  }
}
model {
  //priors

  //One possible regularization
  if(beta_regularization){
    unit_dS[2:n_days] ~ lognormal(log(unit_dS[:n_days-1]), 1 /beta_regularization);
  }
  //This imposes a very wide prior on the proportion of still susceptible people!
  unit_dS[n_days+1] ~ uniform(0, 1);

  omega ~ normal(prior_omega_mean, prior_omega_std);
  dI ~ normal(prior_dI_mean, prior_dI_std);
  dT ~ normal(prior_dT_mean, prior_dT_std);
  phi_inv ~ exponential(5);
  phi_twitter_inv ~ exponential(5);
  twitter_rate ~ exponential(prior_twitter_lambda);

  if (compute_likelihood == 1){
    for (i in 1:n_days) {
      if (use_twitter == 1) {
        symptomaticTweets[i] ~ neg_binomial_2(twitter_rate * state_I[i],
                                              phi_twitter);
      }
    }
    new_weekly_deaths ~ neg_binomial_2(weekly_deaths, phi);
  }
}

generated quantities {
  // Tweets
  vector[n_days] state_tweets;
  for (i in 1:n_days) {
    state_tweets[i] = twitter_rate * state_I[i];
  }

  int pred_deaths[n_days];
  int pred_tweets[n_days];
  for (i in 1:n_days) {
     if (compute_likelihood == 1) {
          pred_deaths[i] = neg_binomial_2_rng(state_D[i], phi);
      }
      if (use_twitter == 1) {
          pred_tweets[i] = neg_binomial_2_rng(twitter_rate *
                                   state_I[i], phi_twitter);
      }
  }
  int pred_weekly_deaths[n_weeks];
  for (i in 1:n_weeks) {
     if (compute_likelihood) {
       pred_weekly_deaths = neg_binomial_2_rng(weekly_deaths, phi);
     }
  }
}
