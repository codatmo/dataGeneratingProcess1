rm(list = ls())
set.seed(43614)
runSim2 = function(runName = runName,
                   nPop = nPop,
                   nDays = nDays,
                   print = PRINT,
                   probInfectPerContact = probInfectPerContact,
                   probTweet = probTweet,
                   probDeath = probDeath) {
  
  dayState = c('i', rep('s', nPop - 1))
  infectionDay = c(1, rep(NA, nPop - 1))
  tweets = rep(0, nPop)
  dead = rep(0,nPop)
  
  meanDaysInfectious = 3
  meanContactsPerDay = 10
  
  #dataM = matrix(nrow = nDays, ncol = 3)
  colNames = c('runName', 'day', 's', 'i', 'r', 'd', 'd_today', 'probInfectPerContact', 'tweets','total_deaths')
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
    df[day,]$probInfectPerContact = probInfectPerContact
    df[day,]$tweets = sum(tweets)
    df[day,]$total_deaths = sum(dead)
    
    #dataM[day,] = c(susceptible, infectious, resolved, dead)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, resolved=%d, died_today=%d, dead=%d, tweets=%d, total_deaths=%d\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$d_today, df[day,]$d, df[day,]$tweets, df[day,]$total_deaths))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'd_today') {
        nextDayState[per] = 'd'
      }
      if (dayState[per] == 'i') {
        if (day - infectionDay[per]  > meanDaysInfectious) {
          nextDayState[per] = 'r'
          if (rbinom(1,1,probDeath) == 1) { #do they die?
            dead[per] = 1;
            nextDayState[per] = 'd_today'
          }
        }
      }
    }
    for (per in 1:nPop) {
      if (dayState[per] == 'i') {
        tweets[per] = rbinom(1, 1, probTweet)
        for (otherPer in sample(1:nPop, meanContactsPerDay)) {
          if (dayState[otherPer] == 's' &&
              rbinom(1, 1, probInfectPerContact) == 1) {
            nextDayState[otherPer] = 'i'
            infectionDay[otherPer] = day
          }
        }
      }
    }
    dayState = nextDayState
  }
  return(df)
}
nWeeks = 10
nDays = nWeeks * 7
nPop = 1e+04
simDf = runSim2(runName = "1/30", nPop = nPop, nDays = nDays, 
                         print = TRUE, 
                         probInfectPerContact = 1/30,
                         probTweet = 1/10,
                         probDeath = 3/100)

deaths = colSums(matrix(simDf$total_deaths, nrow=7))

library("cmdstanr")
library(tidyverse)

# Run from command line: Rscript run.R
# If running from RStudio remember to set the working directory
# >Session>Set Working Directory>To Source File Location

model <- cmdstan_model("stan/tweet_sir_negbin.stan")

# Code modified from 
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

stan_data <- list(n_days = nrow(simDf), y0 = c(nPop -1, 1, 0), t0 = 0, 
                  ts = 1:nDays, N = nPop, cases = simDf$i, 
                  symptomaticTweets =  simDf$tweets, 
                  compute_likelihood = 1,
                  run_twitter = 1,
                  run_SIR = 1)

fit <- model$sample(data=stan_data, output_dir="output", parallel_chains = 4, 
                    chains = 4)

#print(fit$summary(variables = c('pred_cases')))

drawsDf <- fit$draws(variables = c('pred_tweets', 'pred_cases'), format = 'draws_df')

sample_size = 40
sampleDrawsDf <- drawsDf[sample(nrow(drawsDf), sample_size, replace = FALSE),]
colnames(sampleDrawsDf)[1:nDays] <- 1:nDays
#sampleDrawsDf <- cbind(1:nrow(sampleDrawsDf), sampleDrawsDf)
#colnames(sampleDrawsDf)[1] = "draw"
longSampleDrawsDf <- gather(sampleDrawsDf, key=day, value=pred_tweets, 
                            1:all_of(nDays))

longSampleDrawsDf$day <- as.numeric(longSampleDrawsDf$day)
longSampleDrawsDf$y <- as.numeric(longSampleDrawsDf$pred_tweets)

simDf$y = simDf$tweets
simDf$.draw = simDf$runName

# ggplot(data = NULL, aes(day, y, group = .draw)) +
#   geom_line(data = longSampleDrawsDf, alpha = 0.2) +
#   geom_point(data = simDf, color = "orange")


predCasesDf = fit$summary(variables = c('pred_cases'))
predCasesDf$day = 1:nrow(predCasesDf)
predTweetsDf = fit$summary(variables = c('pred_tweets'))
predTweetsDf$day = 1:nrow(predTweetsDf)

ggplot(data = NULL, aes(x = day, y = mean)) +
  geom_ribbon(data = predCasesDf, aes(ymin = q5, ymax = q95, alpha = 0.3)) +
  geom_line(data = predCasesDf) +
  geom_point(data = simDf, aes(y = i), color = "red") +
  geom_ribbon(data = predTweetsDf, aes(ymin = q5, ymax = q95, alpha = 0.3)) +
   geom_line(data = predTweetsDf) + 
   geom_point(data = simDf, aes(y = tweets), color = "blue") +
#  scale_fill_manual(labels = c('infected', 'tweets'))
   labs(y = "mean with 95% interval",
        caption = "predicted shaded, lines with dots for truth") +
   theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
# 
#   