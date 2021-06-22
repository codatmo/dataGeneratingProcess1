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

model <- cmdstan_model("stan/tweet_sird_negbin.stan")

# Code modified from 
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

predicted_category = 'd'
stan_data <- list(n_days = nrow(simDf), y0 = c(nPop - infStartValue, infStartValue, 0, 0), t0 = 0, 
                  ts = 1:nDays, N = nPop, cases = simDf[[predicted_category]], 
                  symptomaticTweets =  simDf$tweets, 
                  compute_likelihood = 1,
                  run_twitter = 0,
                  run_SIR = 1)

fit <- model$sample(data=stan_data, output_dir="output", 
                    parallel_chains = 4, 
                    chains = 4, seed = 4857)

resultsBetaDf = fit$summary(variables = c('beta'))
resultsGammaDf = fit$summary(variables = c('gamma'))
resultsDeathDf = fit$summary(variables = c('deaths'))

cat(paste(sprintf("Simulated vs fit mean values for simulated %s:", 
                   predicted_category),
           sprintf("betaInfRateSim=%.2f, beta=%.2f, sd=%.2f", 
                    betaInfRateSim, resultsBetaDf$mean, resultsBetaDf$sd),
           sprintf("gammaResRateSim=%.2f, gamma=%.2f, sd=%.2f",
                   gammaResRateSim,resultsGammaDf$mean, resultsGammaDf$sd),
           sprintf("deathRateSim=%.3f, deaths=%.3f, sd=%.3f",
                   deathRateSim, resultsDeathDf$mean, resultsDeathDf$sd),
           sep = "\n")
    )

#for a one day delay from d_today -> d, adds extra compartment that is not in 
#the ODE model. Results are:
# Simulated vs fit mean values for simulated d:
# betaInfRateSim=0.30, beta=0.42, sd=0.02
# gammaResRateSim=0.14, gamma=0.26, sd=0.02
# deathRateSim=0.100, deaths=0.004, sd=0.000

# without one day delay d_today = 0, matches ODE compartments:
#Simulated vs fit mean values for simulated d:
#betaInfRateSim=0.30, beta=0.46, sd=0.02
#gammaResRateSim=0.14, gamma=0.30, sd=0.02
#deathRateSim=0.100, deaths=0.004, sd=0.000

predCasesDf = fit$summary(variables = c('pred_cases'))
predCasesDf$day = 1:nrow(predCasesDf)

stan_data2 = stan_data

stan_data2$run_twitter = 1

fitTwitter <- model$sample(data=stan_data2, output_dir="output", 
                           parallel_chains = 4, 
                           chains = 4,
                           seed = 4857)

resultsBetaDf2 = fitTwitter$summary(variables = c('beta'))
resultsGammaDf2 = fitTwitter$summary(variables = c('gamma'))
resultsLambdaTwitter = fitTwitter$summary(variables = c('lambda_twitter'))

cat(paste0(sprintf("Simulated vs fit values for simulated %s:", 
                              predicted_category),
           sprintf("\nbetaInfRateSim=%.2f, beta=%.2f, sd=%.2f", 
                   betaInfRateSim, resultsBetaDf2$mean, resultsBetaDf2$sd),
           sprintf("\ngammaResRateSim=%.2f, gamma=%.2f, sd=%.2f",
                   gammaResRateSim,resultsGammaDf2$mean, resultsGammaDf2$sd)),
           sprintf("\nprobTweetSim=%.2f, lambda_twitter=%.2f, sd=%.2f",
                   probTweetSim, resultsLambdaTwitter$mean, 
                   resultsLambdaTwitter$sd)
)

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
            aes(color = 'predicted cases')) +
  #cases twitter informed
  geom_ribbon(data = predCasesTwitterDf,  aes(ymin = .data[[minQuantileLabel]], 
                                              ymax = .data[[maxQuantileLabel]], 
                                              fill = predCasesRibbon),
              alpha = 0.3, fill = "blue") +
  geom_line(data = predCasesTwitterDf, aes(color = 'predicted cases + twitter')) +
  
  geom_point(data = simDf, aes(y = d, color = 'sim cases'), size = .5, color = "red") +
  
  #twitter
  geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]], 
                                        ymax = .data[[maxQuantileLabel]], 
                                        fill = 'pred tweets'),
              alpha = 0.3, fill = "yellow") +
  geom_line(data = predTweetsDf, aes(color = 'predicted tweets')) + 
  geom_point(data = simDf, aes(y = tweets, color = 'sim tweets'), size = .5, color = 'black') +
  
  #  scale_fill_manual(labels = c('infected', 'tweets'))
  labs(y = "mean with 80% shaded interval and sim data in dots",
       caption = "predicted shaded, lines with dots for simulated truth") +
  theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
