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
  colNames = c('runName', 'day', 's', 'i', 'r', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(colNames)))
  colnames(df) = colNames
  nextDayState = dayState
  for (day in 1:nDays) {
    df[day,] = rep(NA,length(colNames)) #setup data, will cause errors if I miss something
    df[day,]$runName = runName
    df[day,]$s=  length(subset(dayState, dayState == 's'))
    df[day,]$i = length(subset(dayState, dayState == 'i'))
    df[day,]$r = length(subset(dayState, dayState == 'r'))
    df[day,]$day = day
    df[day,]$tweets = sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, resolved=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$tweets, betaInfRate/gammaResRate))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'i') {
        tweets[per] = rbinom(1, 1, probTweet)
        if (rbinom(n = 1, size = 1, prob = gammaResRate) == 1) {  
          nextDayState[per] = 'r'
        }
      }
    }
    for (per in 1:nPop) {
      if (dayState[per] == 'i') {
        for (otherPer in sample(1:nPop, contactPopulationSize)) {
          if (dayState[otherPer] == 's' &&
              rbinom(n = 1, size = 1, prob = betaInfRate/contactPopulationSize) == 1) {
            nextDayState[otherPer] = 'i'
            # lastDayInfected[otherPer] = day + rnorm(1, mean = lambdaResRate, 
            #                                         sd = 2)
          }
        }
      }
    }
    dayState = nextDayState
  }
  return(df)
}




runSim3 = function(runName = runName,
                   nPop = nPop,
                   nDays = nDays,
                   print = PRINT,
                   betaInfRate = betaInfRate,
                   gammaResRate = gammaResRate,
                   probTweet = probTweet,
                   infStartValue = infStartValue) {
  
  dayState = c(rep('i',infStartValue), rep('s', nPop - infStartValue))
  tweets = rep(0, nPop)
  contactPopulationSize = 10
  colNames = c('runName', 'day', 's', 'i', 'r', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(colNames)))
  colnames(df) = colNames
  nextDayState = dayState
  for (day in 1:nDays) {
    df[day,] = rep(NA,length(colNames)) #setup data, will cause errors if I miss something
    df[day,]$runName = runName
    df[day,]$s=  length(subset(dayState, dayState == 's'))
    df[day,]$i = length(subset(dayState, dayState == 'i'))
    df[day,]$r = length(subset(dayState, dayState == 'r'))
    df[day,]$day = day
    df[day,]$tweets = sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, resolved=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$tweets, betaInfRate/gammaResRate))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'i') {
        tweets[per] = rbinom(1, 1, probTweet)
        if (rbinom(n = 1, size = 1, prob = gammaResRate) == 1) {  
          nextDayState[per] = 'r'
        }
      }
    }
    for (per in 1:nPop) {
      if (dayState[per] == 'i') {
        for (otherPer in sample(1:nPop, contactPopulationSize)) {
          if (dayState[otherPer] == 's' &&
              rbinom(n = 1, size = 1, prob = betaInfRate/contactPopulationSize) == 1) {
            nextDayState[otherPer] = 'i'
            # lastDayInfected[otherPer] = day + rnorm(1, mean = lambdaResRate, 
            #                                         sd = 2)
          }
        }
      }
    }
    dayState = nextDayState
  }
  return(df)
}

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
  
  print(sprintf("Beta=%.2f, gamma=%.2f",  meanContactsPerDay * probInfectPerContact,
                meanDaysInfectious))
  
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
        if (day - infectionDay[per]  > rnorm(1, meanDaysInfectious, 1)) {
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
        for (otherPer in sample(1:nPop, rnorm(1, meanContactsPerDay, 3))) {
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
nWeeks = 20
nDays = nWeeks * 7
nPop = 1e+04
# simDf = runSim2(runName = "1/30", nPop = nPop, nDays = nDays, 
#                          print = TRUE, 
#                          probInfectPerContact = 1/30,
#                          probTweet = 1/2,
#                          probDeath = 3/100)


betaInfRateSim = .4
gammaResRateSim = .33
probTweetSim = .5
deathRateSim = .1
infStartValue = 1

simDf = SIRDsim(runName = 'test', nPop = nPop, nDays = nDays, print = TRUE, 
                betaInfRate = betaInfRateSim, gammaResRate = gammaResRateSim, 
                probTweet = probTweetSim, deathRate = deathRateSim, infStartValue = infStartValue)


simDf = runSim3(runName = 'test', nPop = nPop, nDays = nDays, print = TRUE, 
                betaInfRate = betaInfRateSim, gammaResRate = gammaResRateSim, 
                probTweet = probTweetSim, infStartValue = infStartValue)

library("cmdstanr")
library(tidyverse)

# Run from command line: Rscript run.R
# If running from RStudio remember to set the working directory
# >Session>Set Working Directory>To Source File Location

model <- cmdstan_model("stan/tweet_sir_negbin.stan")

# Code modified from 
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

stan_data <- list(n_days = nrow(simDf), y0 = c(nPop - infStartValue, infStartValue, 0), t0 = 0, 
                  ts = 1:nDays, N = nPop, cases = simDf$i, 
                  symptomaticTweets =  simDf$tweets, 
                  compute_likelihood = 1,
                  run_twitter = 0,
                  run_SIR = 1)



fit <- model$sample(data=stan_data, output_dir="output", 
                           parallel_chains = 4, 
                           chains = 4)


resultsBetaDf = fit$summary(variables = c('beta'))
resultsGammaDf = fit$summary(variables = c('gamma'))

cat(paste0(sprintf("Simulated vs fit values: \nbetaInfRateSim=%.2f, beta=%.2f, sd=%.2f", 
                    betaInfRateSim, resultsBetaDf$mean, resultsBetaDf$sd),
             sprintf("\ngammaResRateSim=%.2f, gamma=%.2f, sd=%.2f",
                     gammaResRateSim,resultsGammaDf$mean, resultsGammaDf$sd)))

predCasesDf = fit$summary(variables = c('pred_cases'))
predCasesDf$day = 1:nrow(predCasesDf)

# predCasesSIRDf = fit_no_twitter$summary(variables = c('pred_cases'))
# predCasesSIRDf$day = 1:nrow(predCasesSIRDf)
# predTweetsDf = fit$summary(variables = c('pred_tweets'))
# predTweetsDf$day = 1:nrow(predTweetsDf)

ggplot(data = NULL, aes(x = day, y = mean)) +
geom_ribbon(data = predCasesDf, aes(ymin = q5, ymax = q95, alpha = 0.3)) +
geom_point(data = simDf, aes(y = i), color = "red")

stan_data2 = stan_data

stan_data2$run_twitter = 1
 
fitTwitter <- model$sample(data=stan_data2, output_dir="output", 
                               parallel_chains = 4, 
                               chains = 4)

resultsBetaDf2 = fitTwitter$summary(variables = c('beta'))
resultsGammaDf2 = fitTwitter$summary(variables = c('gamma'))
resultsLambdaTwitter = fitTwitter$summary(variables = c('lambda_twitter'))

cat(paste0(sprintf("Simulated vs fit values: \nbetaInfRateSim=%.2f, beta=%.2f, sd=%.2f", 
                   betaInfRateSim, resultsBetaDf2$mean, resultsBetaDf2$sd),
           sprintf("\ngammaResRateSim=%.2f, gamma=%.2f, sd=%.2f",
                   gammaResRateSim,resultsGammaDf2$mean, resultsGammaDf2$sd)),
           sprintf("\nprobTweetSim=%.2f, lambda_twitter=%.2f, sd=%.2f",
                   probTweetSim, resultsLambdaTwitter$mean, resultsLambdaTwitter$sd)
)

# #print(fit$summary(variables = c('pred_cases')))
# 
# drawsDf <- fit$draws(variables = c('pred_tweets', 'pred_cases'), format = 'draws_df')
# 
# sample_size = 40
# sampleDrawsDf <- drawsDf[sample(nrow(drawsDf), sample_size, replace = FALSE),]
# colnames(sampleDrawsDf)[1:nDays] <- 1:nDays
# #sampleDrawsDf <- cbind(1:nrow(sampleDrawsDf), sampleDrawsDf)
# #colnames(sampleDrawsDf)[1] = "draw"
# longSampleDrawsDf <- gather(sampleDrawsDf, key=day, value=pred_tweets, 
#                             1:all_of(nDays))
# 
# longSampleDrawsDf$day <- as.numeric(longSampleDrawsDf$day)
# longSampleDrawsDf$y <- as.numeric(longSampleDrawsDf$pred_tweets)
# 
# simDf$y = simDf$tweets
# simDf$.draw = simDf$runName
# 
# # ggplot(data = NULL, aes(day, y, group = .draw)) +
# #   geom_line(data = longSampleDrawsDf, alpha = 0.2) +
# #   geom_point(data = simDf, color = "orange")
# 
# 
minQuantile = .2
maxQuantile = .8
minQuantileLabel = '20%'
maxQuantileLabel = '80%'

predCasesDf = fit$summary(variables = c('pred_cases'), mean,
                          ~quantile(.x, probs = c(minQuantile, maxQuantile)))
predCasesDf$day = 1:nrow(predCasesDf)

predCasesTwitterDf = fitTwitter$summary(variables = c('pred_cases'), mean,
                                        ~quantile(.x, 
                                                  probs = c(minQuantile, 
                                                            maxQuantile)))
predCasesTwitterDf$day = 1:nrow(predCasesTwitterDf)


predTweetsDf = fitTwitter$summary(variables = c('pred_tweets'), mean,
                                  ~quantile(.x, 
                                            probs = c(minQuantile, 
                                                      maxQuantile)))
predTweetsDf$day = 1:nrow(predTweetsDf)

predCasesLine = 'predicted cases'
predCasesRibbon = paste0('predicted cases ',minQuantileLabel, '/', maxQuantileLabel)

ggplot(data = NULL, aes(x = day, y = mean)) +
   #cases no twitter info
   geom_ribbon(data = predCasesDf, aes(ymin = .data[[minQuantileLabel]], 
                                       ymax = .data[[maxQuantileLabel]], 
                                       fill = predCasesRibbon), 
               alpha = 0.3, fill = "red") +
   geom_line(data = predCasesDf, 
             aes(color = 'pred cases infection')) +
   #cases twitter informed
   geom_ribbon(data = predCasesTwitterDf,  aes(ymin = .data[[minQuantileLabel]], 
                                               ymax = .data[[maxQuantileLabel]], 
                                               fill = predCasesRibbon),
               alpha = 0.3, fill = "blue") + 
   geom_line(data = predCasesTwitterDf, aes(color = 'pred cases infection + twitter')) +

   geom_point(data = simDf, aes(y = i, color = 'sim cases'), size = .5, color = "black") +

   #twitter
   geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]], 
                                         ymax = .data[[maxQuantileLabel]], 
                                         fill = 'pred tweets'),
               alpha = 0.3, fill = "yellow") +
   geom_line(data = predTweetsDf, aes(color = 'pred tweets')) + 
   geom_point(data = simDf, aes(y = tweets, color = 'sim tweets'), size = .5, color = 'black') +
   
   #  scale_fill_manual(labels = c('infected', 'tweets'))
    labs(y = "mean with 80% shaded interval and sim data in dots",
         caption = "predicted shaded, lines with dots for simulated truth") +
    theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
# # 
# #   