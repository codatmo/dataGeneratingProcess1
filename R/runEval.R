# dependencies
library(tidyverse)
library(cmdstanr)
library(data.table)
library(kableExtra)

source(here::here("R","util.R"))
source(here::here("R","SIRTDsim.R"))
source(here::here("R","sim_configs.R"))
source(here::here("R","modeling_configs.R"))
# dependencies
# setup run_df
run_df <- setup_run_df(seed = 93435, n_pop = 1000, n_days = 7) # in R/util.R
run_df <- sim_brazil_1(run_df) # in R/sim_configs.R
draws_run_df <- sim_draw_params(2, run_df) # in R/sim_configs.R
run_df <- rbind(run_df, draws_run_df) 

run_df <- model_stan_SIRDs(run_df) #in R/modeling_configs.R
run_df$ode_solver <- 'block' # set ODE solver in Stan program across all runs
run_df$compute_likelihood <- 1 # compute likelihood across all runs
run_df$reports <- list(c('graph_sim','graph_ODE', 'graph_tweets', 'graph_d', 'plot','param_recovery'))
# setup run_df
# run models
for (j in 1:nrow(run_df)) {
  fit <- NA
  if (run_df[j,]$model_to_run == 'twitter_sird') {
    stan_data <- 
      list(n_days = run_df[j,]$n_days,
           sDay1 = run_df[j,]$n_pop - 1,
           iDay1 = 1,
           rDay1 = 0,
           dDay1 = 0,            
           NPop = run_df[j,]$n_pop,
           tweets = unlist(run_df[j,]$tweets),
           deaths = unlist(run_df[j,]$d),
           compute_likelihood = run_df[j,]$compute_likelihood,
           run_twitter = run_df[j,]$apply_twitter_data,
           run_block_ODE = ifelse(run_df[j,]$ode_solver == 'block', 1, 0),
           run_rk45_ODE = ifelse(run_df[j,]$ode_solver == 'rk45', 1, 0))
    model <- cmdstan_model(here::here("stan", "tweet_sird.stan"))
    
    fit <- model$sample(data=stan_data,
                        parallel_chains = 4,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        chains = 4,
                        seed = 4857)
    run_df[j,]$fit = list(fit)
  }
  else if (run_df[j,]$model_to_run == 'tweet_sird_negbin_optimized') {
    stan_data_2 <- list(n_days = run_df[j,]$n_days,
                        y0 = c(run_df[j,]$n_pop - run_df[j,]$n_patient_zero, 
                               run_df[j,]$n_patient_zero, 0, 0, 0), # one more zero here
                        t0 = 0,
                        ts = 1:run_df[j,]$n_days,
                        compute_likelihood = run_df[j,]$compute_likelihood,
                        use_twitter = run_df[j,]$apply_twitter_data,
                        death_count = run_df[j,]$d,
                        symptomaticTweets = run_df[j,]$tweets,
                        trapezoidal_solver = 0)
    model2 <- cmdstan_model(here::here("stan", "tweet_sird_negbin_optimized.stan"))
    fit <- model2$sample(data=stan_data_2,
                         parallel_chains = 4,
                         iter_warmup = 1000,
                         iter_sampling = 1000,
                         chains = 4,
                         seed = 4857)
  }
  if (run_df[j,]$model_to_run != 'none') {
    d_tweets_in_interval = countPredictionsInQuantile(fit = fit, 
                                                      run_df = run_df, 
                                                      j = j, print = TRUE)  
    run_df[j,]$d_in_interval = d_tweets_in_interval[1]
    run_df[j,]$tweets_in_interval = d_tweets_in_interval[2]
  }
# run models
  else {
    print(sprintf("no model selected, got:'%s'",run_df[j,]$model_to_run));
  }
# section 6
  plot <- ggplot(data = NULL, aes(x = day, y = count))
  if ('graph_sim' %in% unlist(run_df[j,]$reports)) {
    plot <- plot +  graph_sim_data(data_df = run_df[j,], hide_s = TRUE)
  }
  if ('graph_ODE' %in% unlist(run_df[j,]$reports)) {
    plot <- plot + graph_ODE(data_df = run_df[j,], fit = fit, hide_s = TRUE)
  }
  if ('graph_tweets' %in% unlist(run_df[j,]$reports)) {
    plot <- plot_predictions(plot = plot, prediction_label = 'pred_tweets', 
                             fit = fit, 
                             show_ribbon = TRUE)
  }
  if ('graph_d' %in% unlist(run_df[j,]$reports)) {
    plot <- plot_predictions(plot = plot, prediction_label = 'pred_deaths', 
                             fit = fit, 
                             show_ribbon = TRUE)
  }
  if ('plot' %in% unlist(run_df[j,]$reports)) {
    print(plot)
  }
# section 6
# section 7
  if ('param_recovery' %in% run_df[j,]$reports) {
    cat(param_recovery(data_df = run_df[j,], fit = fit))
  }
# section 7
}
# section 8
summary_cols = c('sim_run_id', 'model_to_run', 'beta_mean', 'gamma', 'death_prob',
                 'tweet_rate', 'days2death', 'ode_solver', 'description',
                 'd_in_interval', 'tweets_in_interval')
print(run_df[,summary_cols])
# section 8
