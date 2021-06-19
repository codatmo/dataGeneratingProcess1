library("cmdstanr")
library(tidyverse)

paramsForUnitTest <- function(runDf, i, seed, nPop, nDays, nDailyContacts,
                             nPatientZero) {
   return(
     sprintf("\n
             set.seed(%d)
             nPop <- %d
             nDays <- %d
             gammaResolvedPerDayRateSim <- %.2f
             tweetRateInfected <- %.2f
             deathRateSim <- %.2f
             meanDaysToDeathFromT <- %.2f
             sdOfBetas <- %.2f
             betaInfectionRateMeanSim <- %.2f
             nPatientZero <- %d
             nDailyContacts <- %d\n",
             seed, nPop, nDays,
             runDf[i,]$gamma,
             runDf[i,]$tweetRate,
             runDf[i,]$deathRate,
             runDf[i,]$daysToDeath, 
             runDf[i,]$sdBeta, 
             runDf[i,]$betaMean,
             runDf[i,]$nPatientZero,
             runDf[i,]$nDailyContacts
             )
   )
}

source(here::here("R","SIRTDsim.R"))
seed = 4614
set.seed(seed)
nPop <- 10000
nWeeks <- 10
nDays <- nWeeks * 7
nSims <- 10

runDf <- data.frame(runId = c(1:nSims))
runDf$gamma <- runif(nSims,0,.5) #c(1/2,1/10)
runDf$tweetRate <-  runif(nSims,0,2) #c(.2,.8)
runDf$deathRate <-  runif(nSims,.01,1) #c(.1,.2)
runDf$daysToDeath <-  runif(nSims,1,20) #c(1,10)
runDf$sdBeta <-  rep(.00001, nSims) #c(.0001,.2)
runDf$betaMean <-  runif(nSims,0,1) #c(.2,.4)
runDf$nPatientZero <-  rep(10,nSims) #c(10,100)
runDf$nDailyContactsInf <- rep(10,nSims)
runDf$betaInfRatePerDay <- list(rep(NA,nDays),rep(NA,nDays))

for (i in 1:nSims) {
  betaForWeek = abs(rnorm(nWeeks, runDf[i,]$betaMean, runDf[i,]$sdBeta))
  runDf[i,]$betaInfRatePerDay = list(rep(betaForWeek, times = rep(7,nWeeks)))
  summary = sprintf("%d, B=%.2f, Bsd=%.2f, G=%.2f, D=%.2f, #D =%.2f, T=%.2f",runDf[i,]$runId, 
                    runDf[i,]$betaMean, runDf[i,]$sdBeta, runDf[i,]$gamma,
                    runDf[i,]$deathRate, runDf[i,]$daysToDeath, 
                    runDf[i,]$tweetRate)
  
  if (FALSE) {
    cat(paramsForUnitTest(runDf, i, seed, nPop, nDays))
  }
  
  simDf = SirtdVaryBeta(runName = runDf[i,]$runId, 
                        nPop = nPop, 
                        nDays = nDays, 
                        print = TRUE,
                        betaDailyInfectionRates = unlist(runDf[i,]$betaInfRatePerDay),
                        gammaResolvedPerDayRate = runDf[i,]$gamma,
                        tweetRateInfected = runDf[i,]$tweetRate,
                        meanDaysToDeathFromT = runDf[i,]$daysToDeath,
                        nPatientZero = runDf[i,]$nPatientZero,
                        nDailyContacts = runDf[i,]$nDailyContactsInf,
                        deathRate = runDf[i,]$deathRate)

  simData <- matrix(data = c(simDf$s, simDf$i, simDf$r, simDf$d,
                          simDf$tweets),
                 nrow = nrow(simDf), ncol = 5)
	compartment = 4
  tweetSourceIndex <- 2
  stan_data <- list(n_days = nrow(simDf),
                  y0 = c(nPop - runDf[i,]$nPatientZero, runDf[i,]$nPatientZero, 0, 0),
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
                  run_twitter = 1,
                  run_SIR = 1,
                  check_ODE = 0)

		  model <- cmdstan_model(here::here("stan", "tweet_sird_negbin.stan"))

      fit <- model$sample(data=stan_data,
                    parallel_chains = 4,
                    iter_warmup = 1000,
                    iter_sampling = 1000,
                    chains = 4,
                    seed = 4857)
      
      minQuantile = .2
      maxQuantile = .8
      minQuantileLabel = '20%'
      maxQuantileLabel = '80%'
      
      predCasesDf = fit$summary(variables = c('pred_cases'), mean,
                                ~quantile(.x, probs = c(minQuantile, maxQuantile)))
      predCasesDf$day = 1:nrow(predCasesDf)
      
      predTweetsDf = fit$summary(variables = c('pred_tweets'), mean,
                                 ~quantile(.x, probs = c(minQuantile, maxQuantile)))
      predTweetsDf$day = 1:nrow(predTweetsDf)
      
      
      predCasesLine = 'predicted cases'
      predCasesRibbon = paste0('predicted cases ',minQuantileLabel, '/', maxQuantileLabel)
      
      plot = ggplot(data = NULL, aes(x = day, y = mean)) +
        #cases no twitter info
        geom_ribbon(data = predCasesDf, aes(ymin = .data[[minQuantileLabel]], 
                                            ymax = .data[[maxQuantileLabel]], 
                                            fill = predCasesRibbon), 
                    alpha = 0.3, fill = "red") +
        geom_line(data = predCasesDf, 
                  aes(color = 'pred cases death')) +
        #cases twitter informed
        geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]], 
                                              ymax = .data[[maxQuantileLabel]], 
                                              fill = predCasesRibbon),
                    alpha = 0.3, fill = "blue") + 
        geom_point(data = simDf, aes(y = d, color = 'sim cases'), size = .5, color = "black") +
        
        #twitter
        geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]], 
                                              ymax = .data[[maxQuantileLabel]], 
                                              fill = 'pred tweets'),
                    alpha = 0.3, fill = "yellow") +
        geom_line(data = predTweetsDf, aes(color = 'pred tweets')) + 
        geom_point(data = simDf, aes(y = tweets, color = 'sim tweets'), size = .5, color = 'black') +
        
        #  scale_fill_manual(labels = c('infected', 'tweets'))
        labs(y = "mean with 80% shaded interval and sim data in dots",
             caption = summary) +
        theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
      print(plot)
        
}
    
