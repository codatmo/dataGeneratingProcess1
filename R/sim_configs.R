source(here::here("R","util.R"))
source(here::here("R","SIRTDsim.R"))

# Each run is a row in runDf, the rows contain all information necessary to run
# an experiment. The data in the columns can vary based on whether a sim is being
# run, real data etc. 
# Brazil data
# Simulated data
# 
# section 1
#Globals

#' Returns a value +- the variation param around a mean value.
#' @param mean 
#' @param variation
#' @return 
draw_unif_prob <- function(mean, variation) {
  if (mean > 1 || mean < 0) {
    stop(sprintf("mean out of range for probability in draw_unif_prob %d",
               mean))
  }
  return (round(runif(n = 1, 
          min = max(mean - variation, 0), 
          max = min(mean + variation, 1)),3))
}

#' Takes source_df with values for simulation parameters and draws n_sims 
#' times from uniform values +- .2 for Beta, Lambda, .05 for death_prob,
#' and so on, look to source for details in 'R/data_config.R'
#' @param n_sims number of draws for params
#' @param source_df run_df with point estimates that will center draws
#' @return run_df with n_sims number of draws, source_df not included
sim_draw_params <- function(n_sims, source_df) {
  last_sim_id <- source_df$sim_run_id  
  return_df <- data.frame()
  for (j in 1:n_sims) {
    drawn_df <- copy(source_df)
    last_sim_id <- last_sim_id + 1
    drawn_df$sim_run_id <- last_sim_id
    drawn_df$beta_mean <- draw_unif_prob(mean = source_df$beta_mean,
                                         variation = 0.2)
    drawn_df$beta_daily_rate <- list(rep(drawn_df$beta_mean, 
                                         drawn_df$n_days))
    drawn_df$gamma <- draw_unif_prob(mean = source_df$gamma,
                                     variation = 0.2)
    drawn_df$tweet_rate <- draw_unif_prob(mean = source_df$tweet_rate,
                                          variation = .25)
    drawn_df$death_prob <- draw_unif_prob(mean = source_df$death_prob,
                                          variation = .02)
    sim_df <- sirtd_vary_beta_exact(seed = drawn_df$seed,
                                    n_pop = drawn_df$n_pop, 
                                    n_days = drawn_df$n_days,
                                    print = FALSE,
                                    beta_daily_inf_rates = 
                                      unlist(drawn_df$beta_daily_rate),
                                    gamma_res_per_day_rate = drawn_df$gamma,
                                    tweet_rate_infected = drawn_df$tweet_rate,
                                    mean_days_to_death_from_t = drawn_df$days2death,
                                    n_patient_zero = drawn_df$n_patient_zero,
                                    death_prob = drawn_df$death_prob)
    drawn_df$s <- list(sim_df$s)
    drawn_df$i <- list(sim_df$i)
    drawn_df$r <- list(sim_df$r)
    drawn_df$t <- list(sim_df$t)
    drawn_df$d <- list(sim_df$d)
    drawn_df$tweets <- list(sim_df$tweets)
    return_df <- rbind(return_df, drawn_df)
  }
  return(return_df)
}

#' Mirrors Brazil's observerd deaths/tweets roughly with simulation and given 
#' params, no model configuration/reporting configured
#' @return dataframe for use in runEval.R
sim_brazil_1 <- function (source_df) {
  run_df <- copy(source_df)
  run_df$beta_mean <- .30
  run_df$beta_daily_rate <- list(rep(run_df$beta_mean, 
                                          run_df$n_days))
  run_df$gamma <- .255
  run_df$death_prob <- .01
  run_df$tweet_rate <- .1
  run_df$days2death <- 10
  run_df$n_patient_zero <- 10
  run_df$description <- "Brazil 1 approximation"
  run_df$sim_run_id <- 1
  sim_df <- sirtd_vary_beta_exact(seed = run_df$seed,
                            n_pop = run_df$n_pop, 
                            n_days = run_df$n_days,
                            print = FALSE,
                            beta_daily_inf_rates = 
                              unlist(run_df$beta_daily_rate),
                            gamma_res_per_day_rate = run_df$gamma,
                            tweet_rate_infected = run_df$tweet_rate,
                            mean_days_to_death_from_t = run_df$days2death,
                            n_patient_zero = run_df$n_patient_zero,
                            death_prob = run_df$death_prob)
  run_df$s <- list(sim_df$s)
  run_df$i <- list(sim_df$i)
  run_df$r <- list(sim_df$r)
  run_df$t <- list(sim_df$t)
  run_df$d <- list(sim_df$d)
  run_df$tweets <- list(sim_df$tweets)
  return(run_df)
}
#214,110,287
#run_df <- setup_run_df(seed = 123, n_pop = 214110287, n_days = 365)
#sim_brazil_1(run_df)
