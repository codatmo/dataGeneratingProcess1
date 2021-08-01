#todo
# serialize run_df?? json?
# run Brazil data with baseline to get params 
#
# dependencies
library(tidyverse)
library(cmdstanr)
library(data.table)
library(kableExtra)

source(here::here("R","util.R"))
source(here::here("R","SIRTDsim.R"))
source(here::here("R","sim_configs.R"))
source(here::here("R", "data_configs.R"))
source(here::here("R","modeling_configs.R"))
# dependencies
# setup run_df
#run_df <- setup_run_df(seed = 93435, n_pop = 214110287, n_days = 300) # in R/util.R
run_df <- setup_run_df(seed = 93435, n_pop = 214110287, n_days = 300) # in R/util.R
brazil_sim_df <- sim_brazil_1(run_df) # in R/sim_configs.R
brazil_actual_df <- data_brazil_1(run_df)
#draws_run_df <- sim_draw_params(2, run_df) # in R/sim_configs.R
run_df <- rbind(brazil_sim_df, brazil_actual_df)

#run_df <- model_stan_baseline(run_df) #in R/modeling_configs.R
#run_df$compute_likelihood <- 1 # compute likelihood across all runs
#run_df$reports <- list(c('graph_sim', 'graph_tweets', 'graph_d', 'plot'))
   #list(c('graph_sim','graph_ODE', 'graph_tweets', 'graph_d', 'plot','param_recovery'))
# setup run_df
# run models
j <- 0
while (j < nrow(run_df)) {
  j <- j + 1
  fit <- NA
  if (run_df[j,]$model_to_run == 'baseline') {
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
           run_rk45_ODE = ifelse(run_df[j,]$ode_solver == 'rk45', 1, 0),
           scale = 1,
           center = 1)
    model <- cmdstan_model(here::here("stan", "baseline.stan"))

    fit <- model$sample(data=stan_data,
                        parallel_chains = 4,
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        chains = 4,
                        seed = 4857)
    run_df[j,]$fit = list(fit)
  }
  else if (run_df[j,]$model_to_run == 'UNINOVE_Brazil') {
    stan_data_2 <- list(n_days = run_df[j,]$n_days,
                        y0 = c(run_df[j,]$n_pop - run_df[j,]$n_patient_zero,
                               run_df[j,]$n_patient_zero, 0, 0, 0), # one more zero here
                        t0 = 0,
                        ts = 1:run_df[j,]$n_days,
                        compute_likelihood = run_df[j,]$compute_likelihood,
                        use_twitter = run_df[j,]$apply_twitter_data,
                        death_count = unlist(run_df[j,]$d),
                        symptomaticTweets = unlist(run_df[j,]$tweets),
                        prior_beta_mean = 0,
                        prior_beta_std = 10,
                        prior_omega_mean = 0,
                        prior_omega_std = 10,
                        prior_dI_mean = 0,
                        prior_dI_std = 10,
                        prior_dT_mean = 0,
                        prior_dT_std = 10,
                        prior_twitter_lambda = 1
                        )
    model2 <- cmdstan_model(here::here("stan", "tweet_sirtd_negbin_ODE.stan"))
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
  if ('graph_data' %in% unlist(run_df[j,]$reports)) {
    plot <- graph_real_data(data_df = run_df[j,], plot = plot)
  }
  if ('graph_sim' %in% unlist(run_df[j,]$reports)) {
    plot <- graph_sim_data(data_df = run_df[j,], hide_s = TRUE, plot = plot)
  }
  if ('graph_ODE' %in% unlist(run_df[j,]$reports)) {
    plot <- graph_ODE(data_df = run_df[j,], fit = fit, hide_s = TRUE,
                             plot = plot)
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
  plot = plot + theme(legend.position = "none")
  if (length(unlist(run_df[j,]$reports)) > 0) {
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
                 'tweet_rate', 'days2death', 'description',
                 'd_in_interval', 'tweets_in_interval','n_days')
print(run_df[,summary_cols])
# section 8
