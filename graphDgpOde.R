rm(list = ls())
set.seed(43614)

SIRDsim = function(runName = runName,
                   nPop = nPop,
                   nDays = nDays,
                   print = PRINT,
                   betaInfRate = betaInfRate,
                   gammaResRate = gammaResRate,
                   deathRate = deathRate,
                   probTweet = probTweet,
                   infStartValue = infStartValue) {
  
  dayState = c(rep('i',infStartValue), rep('s', nPop - infStartValue))
  tweets = rep(0, nPop)
  contactPopulationSize = 10
  colNames = c('runName', 'day', 's', 'i', 'r', 'd_today', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(colNames)))
  colnames(df) = colNames
  nextDayState = dayState
  for (day in 1:nDays) {
    df[day,] = rep(NA,length(colNames)) #setup data, will cause errors if I miss something
    df[day,]$runName = runName
    df[day,]$s=  length(subset(dayState, dayState == 's'))
    df[day,]$i = length(subset(dayState, dayState == 'i'))
    df[day,]$r = length(subset(dayState, dayState == 'r'))
    df[day,]$d_today = length(subset(dayState, dayState == 'd_today'))
    df[day,]$d = length(subset(dayState, dayState == 'd'))
    df[day,]$day = day
    df[day,]$tweets = sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, resolved=%d, dead today=%d, 
          dead=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$d_today, 
          df[day,]$d, df[day,]$tweets, betaInfRate/gammaResRate))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'i') {
        tweets[per] = rbinom(1, 1, probTweet)
        if (rbinom(n = 1, size = 1, prob = gammaResRate) == 1) {  
          if (rbinom(n = 1, size = 1, prob = deathRate) == 1) {
            nextDayState[per] = 'd' #not doing daily *1
          }
          else {
            nextDayState[per] = 'r'
          }
        }
      }
      if (dayState[per] == 'd_today') { #irrelvant if not doing daily, see *1
        nextDayState[per] = 'd'
      }
    }
    for (per in 1:nPop) {
      if (dayState[per] == 'i') {
        for (otherPer in sample(1:nPop, contactPopulationSize)) {
          if (dayState[otherPer] == 's' &&
              rbinom(n = 1, size = 1, prob = betaInfRate/contactPopulationSize) == 1) {
            nextDayState[otherPer] = 'i'
          }
        }
      }
    }
    dayState = nextDayState
  }
  return(df)
}

nPop = 10000
nWeeks = 10
nDays = nWeeks * 7
betaInfRateSim = .3
gammaResRateSim = 1/7
probTweetSim = .5
deathRateSim = .1
infStartValue = 10

simDf = SIRDsim(runName = 'test', nPop = nPop, nDays = nDays, print = TRUE, 
                betaInfRate = betaInfRateSim, gammaResRate = gammaResRateSim, 
                probTweet = probTweetSim, deathRate = deathRateSim, 
                infStartValue = infStartValue)


library("cmdstanr")
library(tidyverse)

# Run from command line: Rscript run.R
# If running from RStudio remember to set the working directory
# >Session>Set Working Directory>To Source File Location

model <- cmdstan_model("/home/breck/git/codatmo/dataGeneratingProcess1/stan/tweet_sird_negbin.stan")

# Code modified from 
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

simData = matrix(data = c(simDf$s, simDf$i, simDf$r, simDf$d, 
                          simDf$tweets), 
                 nrow = nrow(simDf), ncol = 5)

compartment = 4
tweetSourceIndex = 2
stan_data <- list(n_days = nrow(simDf), 
                  y0 = c(nPop - infStartValue, infStartValue, 0, 0), 
                  t0 = 0, 
                  ts = 1:nDays, 
                  N = nPop, 
                  n_compartments = ncol(simData) - 1,
                  nDataCols = ncol(simData),
                  compartmentDays = simData, 
                  compartment = compartment, 
                  tweetIndex = 5,
                  tweetSourceIndex = tweetSourceIndex,
                  compute_likelihood = 1,
                  run_twitter = 0,
                  run_SIR = 1)

fit <- model$sample(data=stan_data, 
                    parallel_chains = 1,
                    iter_warmup = 1000,
                    iter_sampling = 100,
                    chains = 1,
                    seed = 4857)

resultsOdeDf = fit$summary(variables = c('ode_states'))
compartmentNames = c('s','i','r','d');
compartmentName = compartmentNames[compartment]
tweetSourceName = compartmentNames[tweetSourceIndex]
odeDf = data.frame(matrix(resultsOdeDf$median, nrow = nrow(simDf), 
       ncol = length(compartmentNames)))
colnames(odeDf) = compartmentNames
odeDf$day = 1:nrow(simDf)

odeLongDf = gather(data = odeDf, key = "compartmentODE", value = "mean",
                   all_of(compartmentNames))

simLongDf = gather(data = simDf, key = "compartmentSim", value = "count",
                   all_of(c('tweets', compartmentNames)))

ggplot(data = NULL, aes(x = day, y = mean)) +
  geom_line(data = odeLongDf, aes(color = compartmentODE)) +
  geom_point(data = simLongDf, aes(y = count, color = compartmentSim), size = .5) +
  labs(y = "median with sim data in dots",
     caption = paste0("predicted shaded, lines with dots for simulated truth, compartment is = ",
     compartmentName), " tweets from ", tweetSourceName) +
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
