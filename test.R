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
  deaths = rep(0,nPop)
  
  meanDaysInfectious = 3
  meanContactsPerDay = 10
  
  #dataM = matrix(nrow = nDays, ncol = 3)
  colNames = c('runName', 'day', 's', 'i', 'r', 'probInfectPerContact', 'tweets','deaths')
  df = data.frame(matrix(nrow = 0, ncol=length(colNames)))
  colnames(df) = colNames
  nextDayState = dayState
  for (d in 1:nDays) {
    df[d,] = rep(NA,length(colNames)) #setup data, will cause errors if I miss something
    df[d,]$runName = runName
    df[d,]$s=  length(subset(dayState, dayState == 's'))
    df[d,]$i = length(subset(dayState, dayState == 'i'))
    df[d,]$r = length(subset(dayState, dayState == 'r'))
    df[d,]$day = d
    df[d,]$probInfectPerContact = probInfectPerContact
    df[d,]$tweets = sum(tweets)
    df[d,]$deaths = sum(deaths)
    
    #dataM[d,] = c(susceptible, infectious, resolved)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, resolved=%d, tweets=%d, deaths=%d\n",
          df[d,]$day, df[d,]$s, df[d,]$i, df[d,]$r, df[d,]$tweets, df[d,]$deaths))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'i') {
        if (d - infectionDay[per]  > meanDaysInfectious) {
          nextDayState[per] = 'r'
          if (rbinom(1,1,probDeath) == 1) { #do they die?
            deaths[per] = 1;
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
            infectionDay[otherPer] = d
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
deathsTweetsDf = runSim2(runName = "1/30", nPop = nPop, nDays = nDays, print = TRUE, 
                         probInfectPerContact = 1/30,
                         probTweet = 1/10,
                         probDeath = 3/100)

deaths = colSums(matrix(deathsTweetsDf$deaths, nrow=7))

library("cmdstanr")
library(tidyverse)

# Run from command line: Rscript run.R
# If running from RStudio remember to set the working directory
# >Session>Set Working Directory>To Source File Location

model <- cmdstan_model("stan/tweet_sir_negbin.stan")

stan_data <- list(n_days = nrow(deathsTweetsDf), y0 = c(nPop -1, 1, 0), t0 = 0, 
                  ts = 1:nDays, N = nPop, cases = deathsTweetsDf$i, 
                  symptomaticTweets =  deathsTweetsDf$tweets, 
                  compute_likelihood = 1,
                  run_twitter = 1,
                  run_SIR = 1)

fit <- model$sample(data=stan_data, output_dir="output", parallel_chains = 4, 
                    chains = 4)
print(fit$summary(variables = c('pred_cases')))

drawsDf <- fit$draws(variables = c('pred_cases'), format = 'draws_df')

sample_size = 40
sampleDrawsDf <- drawsDf[sample(nrow(drawsDf), sample_size, replace = FALSE),]
colnames(sampleDrawsDf)[1:nDays] <- 1:nDays
#sampleDrawsDf <- cbind(1:nrow(sampleDrawsDf), sampleDrawsDf)
#colnames(sampleDrawsDf)[1] = "draw"
longSampleDrawsDf <- gather(sampleDrawsDf, key=day, value=pred_infected, 
                            1:all_of(nDays))

longSampleDrawsDf$day <- as.numeric(longSampleDrawsDf$day)
longSampleDrawsDf$y <- as.numeric(longSampleDrawsDf$pred_infected)

deathsTweetsDf$y = deathsTweetsDf$i
deathsTweetsDf$.draw = deathsTweetsDf$runName


ggplot(data = NULL, aes(day, y, group = .draw)) +
  geom_line(data = longSampleDrawsDf, alpha = 0.2) +
  geom_point(data = deathsTweetsDf, color = "red")


  
  
  