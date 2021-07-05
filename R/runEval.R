library(tidyverse)
library(cmdstanr)
library(data.table)
library(kableExtra)
rm(list = ls())
loadBrazil = FALSE
loadSim = TRUE
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
seed = 93435
n_pop <- 20000
n_days <- 70

# Set up our template, all runs need this info
# Each row of dataframe is a separate run
template_df <- data.frame(sim_run_id = c(NA))
template_df$description <- NA
template_df$seed <- seed

# any simulation/model run
template_df$n_pop <- n_pop
template_df$n_days <- n_days

# any model setup
template_df$model_to_run <- 'none'
template_df$compute_likelihood <- NA

#setup data columns 
template_df$s <- list(c())
template_df$i <- list(c())
template_df$r <- list(c())
template_df$t <- list(c())
template_df$d <- list(c())
template_df$tweets <- list(c())

# setup prediction columns
template_df$d_in_interval <- NA_integer_
template_df$tweets_in_interval <- NA_integer_

# section 1
# section 2
# setup for current simulations

set.seed(template_df$seed)
n_runs <- 2
SIRTD_sim_df = copy(template_df)
if (n_runs > 1) {
  for (j in 2:n_runs) {
    SIRTD_sim_df <- rbind(SIRTD_sim_df,copy(template_df))
  }
}

SIRTD_sim_df$beta_mean <- runif(n_runs, 0.2, 0.3)
SIRTD_sim_df$beta_daily_rate <- list(rep(SIRTD_sim_df$beta_mean, 
                                    SIRTD_sim_df$n_days))
SIRTD_sim_df$gamma <- runif(n_runs, 0.2, 0.3)
SIRTD_sim_df$death_prob <- runif(n_runs, 0.01, 0.03)
SIRTD_sim_df$tweet_rate <- runif(n_runs, 0.1, 1.5)
SIRTD_sim_df$days2death <- 20
SIRTD_sim_df$n_patient_zero <- 20
SIRTD_sim_df$description <- sprintf("SIRTD sim: %d runs", n_runs)
SIRTD_sim_df$sim_run_id <- 1:n_runs
# section 2
# section 3
for (j in 1:n_runs) {
  sim_df <- sirtd_vary_beta(seed = SIRTD_sim_df[j,]$seed,
                         n_pop = SIRTD_sim_df[j,]$n_pop, 
                         n_days = SIRTD_sim_df[j,]$n_days,
                         print = FALSE,
                         beta_daily_inf_rates = 
                                       unlist(SIRTD_sim_df[j,]$beta_daily_rate),
                         gamma_res_per_day_rate = SIRTD_sim_df[j,]$gamma,
                         tweet_rate_infected = SIRTD_sim_df[j,]$tweet_rate,
                         mean_days_to_death_from_t = SIRTD_sim_df[j,]$days2death,
                         n_patient_zero = SIRTD_sim_df[j,]$n_patient_zero,
                         death_prob = SIRTD_sim_df[j,]$death_prob)
  SIRTD_sim_df[j,]$s <- list(sim_df$s)
  SIRTD_sim_df[j,]$i <- list(sim_df$i)
  SIRTD_sim_df[j,]$r <- list(sim_df$r)
  SIRTD_sim_df[j,]$t <- list(sim_df$t)
  SIRTD_sim_df[j,]$d <- list(sim_df$d)
  SIRTD_sim_df[j,]$tweets <- list(sim_df$tweets)
}
# section 3
# section 4
# setup for current model
model_no_tweets_df <- copy(SIRTD_sim_df)
model_no_tweets_df$apply_twitter_data <- 0
model_no_tweets_df$model_to_run <- 'twitter_sird'
model_no_tweets_df$description <- paste0(model_no_tweets_df$description,
                                         " no tweets")

model_tweets_df <- copy(SIRTD_sim_df)
model_tweets_df$apply_twitter_data <- 1
model_tweets_df$model_to_run <- 'twitter_sird'
model_tweets_df$description <- paste0(model_tweets_df$description,
                                      " tweets")
run_df <- rbind(model_no_tweets_df, model_tweets_df)
run_df$ode_solver <- 'block'
run_df$compute_likelihood <- 1
# section 4


# if (loadBrazil) {
#   tweets = read.csv(here::here("data","tweet_count.csv"))
#   colnames(tweets) = c('dateT','predicted')
#   brazilDeaths = readRDS(here::here("data","brazil_nation.rds"))
# 
#   dataDf = data.frame(date = brazilDeaths[1:354,]$date, 
#                     d = brazilDeaths[1:354,]$last_available_deaths)
#   tweetsPadded = rbind(data.frame(dateT = rep(NA,106), 
#                                 predicted = rep(0,106)), 
#                      tweets)
#   dataDf = cbind(dataDf,tweetsPadded)
#   colnames(dataDf) = c('date','d','dateT','tweets')
#   dataDf = dataDf %>% mutate(perc_d = d/nPop) %>% 
#     mutate(perc_t = tweets/nPop)
#   dataDf = dataDf[1:nDays,] #doing 2020 only
# }

# section 5
for (j in 1:nrow(run_df)) {
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
    else {
      print(sprintf("no model selected, got:'%s'",run_df[j,]$model_to_run));
    }
    d_tweets_in_interval = countPredictionsInQuantile(fit = fit, run_df = run_df, j = j, print = TRUE)  
    run_df[j,]$d_in_interval = d_tweets_in_interval[1]
    run_df[j,]$tweets_in_interval = d_tweets_in_interval[2]
}
# section 5

# section 6
  summary_cols = c('sim_run_id', 'model_to_run', 'beta_mean', 'gamma', 'death_prob',
                   'tweet_rate', 'days2death', 'ode_solver', 'description',
                   'd_in_interval', 'tweets_in_interval')
  print(run_df[,summary_cols])
# section 6

graphSimData <- function() {
  if (run_df[j,]$model_to_run == 'none') {
    sim_df = data.frame(day = 1:run_df[j,]$n_days, 
                        tweets = unlist(run_df[j,]$tweets), 
                        s = unlist(run_df[j,]$s),
                        i = unlist(run_df[j,]$i),
                        r = unlist(run_df[j,]$r),
                        t = unlist(run_df[j,]$t),
                        d = unlist(run_df[j,]$d))
    sim_long_df = gather(data = sim_df, key = "compartment_sim", value = "count",
                         all_of(c('tweets', compartmentNames)))
    
    plot = ggplot(data = NULL, aes(x = day, y = count)) +
      geom_point(data = sim_long_df, aes(y = count, 
                                         color = compartment_sim)) +
      labs(y = "sim data in dots",
           caption = paste0("dots for simulated truth") +
             theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0))))
    
    print(plot)
}

  
  if (FALSE) {
    
    
    
    
    plot <- ggplot(data = NULL, aes(x = day, y = mean)) +
      # geom_ribbon(data = predCasesDf, aes(ymin = .data[[minQuantileLabel]],
      #                                     ymax = .data[[maxQuantileLabel]],
      #                                     fill = predCasesRibbon),
      #             alpha = 0.3, fill = "red") +
      geom_line(data = predCasesDf,
                aes(color = 'pred cases death')) +
      geom_point(data = dataDf, 
                 aes(y = perc_d, color = 'sim cases'), size = .5, color = "black")
    
    #   
    #           #cases twitter informed
    # geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]],
    #                                       ymax = .data[[maxQuantileLabel]],
    #                                       fill = predCasesRibbon),
    #             alpha = 0.3, fill = "blue") +
    # 
    # 
    # #twitter
    # geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]],
    #                                       ymax = .data[[maxQuantileLabel]],
    #                                       fill = 'pred tweets'),
    #             alpha = 0.3, fill = "yellow") +
    # geom_line(data = predTweetsDf, aes(color = 'pred tweets')) +
    # geom_point(data = dataDf, aes(y = tweets, color = 'sim tweets'), size = .5, color = 'black') +
    # 
    # #  scale_fill_manual(labels = c('infected', 'tweets'))
    # labs(y = "mean with 80% shaded interval and sim data in dots",
    #      caption = summary) +
    # theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
    print(plot)
  }
  if (loadSim && FALSE) {
    resultsOdeDf = fit$summary(variables = c('ode_states'))
    
    tweetSourceName = compartmentNames[2]
    odeDf = data.frame(matrix(resultsOdeDf$median, nrow = nrow(simDf),
                              ncol = length(compartmentNames)))
    colnames(odeDf) = compartmentNames
    odeDf$day = 1:nrow(simDf)
    
    odeLongDf = gather(data = odeDf, key = "compartmentODE", value = "mean",
                       all_of(compartmentNames))
    
    simLongDf = gather(data = simDf, key = "compartmentSim", value = "count",
                       all_of(c('tweets', compartmentNames)))
    
    plot = ggplot(data = NULL, aes(x = day, y = nPop * mean)) +
      geom_line(data = odeLongDf, aes(color = compartmentODE)) +
      geom_line(data = simLongDf, linetype = 'dotted', aes(y = count * reduction, 
                                                           color = compartmentSim)) +
      labs(y = "median with sim data in dots",
           caption = paste0("predicted shaded, lines with dots for simulated truth, compartment is = ",
                            compartmentName), " tweets from ", tweetSourceName) +
      theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
    
    print(plot)
  }
  
  if (loadSim && FALSE) {
    recovPars = 
      fit$summary(variables = c('gamma', 'beta', 'deathRate', 'lambda_twitter'), 
                  mean, sd)
    cat(sprintf("\nBeta sim = %.4f vs recovered %.4f, sd=%.4f",
                runDf[i,]$betaMean, 
                recovPars[recovPars$variable == 'beta',]$mean,
                recovPars[recovPars$variable == 'beta',]$sd))
    
    cat(sprintf("\nGamma sim = %.4f vs recovered %.4f, sd=%.4f",
                runDf[i,]$gamma, 
                recovPars[recovPars$variable == 'gamma',]$mean,
                recovPars[recovPars$variable == 'gamma',]$sd))
    
    cat(sprintf("\nDeaths sim = %.4f vs recovered %.4f, sd=%.4f",
                runDf[i,]$deathRate, 
                recovPars[recovPars$variable == 'deaths',]$mean,
                recovPars[recovPars$variable == 'deaths',]$sd))
    cat(sprintf("\nLambda Twitter sim = %.4f vs recovered %.4f, sd=%.4f",
                runDf[i,]$tweetRate, 
                recovPars[recovPars$variable == 'lambda_twitter',]$mean,
                recovPars[recovPars$variable == 'lambda_twitter',]$sd))
  }
#  saveRDS(runDf,here::here("R",sprintf("%d_of_%devalBrazil622.rds",i,nrow(runDf))))
}

