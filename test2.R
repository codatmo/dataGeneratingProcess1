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
  colNames = c('runName', 'day', 's', 'i', 'r', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(colNames)))
  colnames(df) = colNames
  nextDayState = dayState
  for (day in 1:nDays) {
    df[day,] = rep(NA,length(colNames)) #setup data, will cause errors if I miss something
    df[day,]$runName = runName
    df[day,]$s=  length(subset(dayState, dayState == 's'))
    df[day,]$i = length(subset(dayState, dayState == 'i'))
    df[day,]$r = length(subset(dayState, dayState == 'r'))
    df[day,]$d = length(subset(dayState, dayState == 'd'))
    df[day,]$day = day
    df[day,]$tweets = sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, resolved=%d, dead=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$d, 
          df[day,]$tweets, betaInfRate/gammaResRate))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'i') {
        tweets[per] = rbinom(1, 1, probTweet)
        if (rbinom(n = 1, size = 1, prob = gammaResRate) == 1) {  
          if (rbinom(n = 1, size = 1, prob = deathRate) == 1) {
            nextDayState[per] = 'd'
          }
          else {
            nextDayState[per] = 'r'
          }
        }
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

nPop = 1000
nWeeks = 10
nDays = nWeeks * 7
betaInfRateSim = .3
gammaResRateSim = .25
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

model <- cmdstan_model("stan/tweet_sird_negbin.stan")

# Code modified from 
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

stan_data <- list(n_days = nrow(simDf), y0 = c(nPop - infStartValue, infStartValue, 0, 0), t0 = 0, 
                  ts = 1:nDays, N = nPop, cases = simDf$i, 
                  symptomaticTweets =  simDf$tweets, 
                  compute_likelihood = 1,
                  run_twitter = 1,
                  run_SIR = 1)



fit <- model$sample(data=stan_data, output_dir="output", 
                    parallel_chains = 4, 
                    chains = 4)


