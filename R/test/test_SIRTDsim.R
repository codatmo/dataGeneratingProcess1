library(testthat)

source(here::here('R','SIRTDsim.R'))

test_that("Tweet test", {
  nPop <- 210147125
  nDays <- 311
  nSims <- 1
  nConfigs <- 1
  runDf <- data.frame(runId = NA)
  runDf$description <- NA
  runDf$tweetsEvalInInterval <- NA
  runDf$casesEvalInInterval <- NA
  runDf$odeSolver <- NA
  runDf$runTwitter <- NA
  runDf$nPop <- nPop
  runDf$nDays <- nDays
  runDf$seed <- 23454
  runDf$modelToRun <- NA
  runDf$compute_likelihood <- 1
  
  # Simulator params
  # replicate rows for nSims
  nSims <- 1
  simRunDf = runDf
  simRunDf$betaMean <-  .3 #runif(nSims, 0, 1) #c(.2,.4)
  simRunDf$gamma <- .29 #runif(nSims, 0, 1) #c(1/2,1/10)
  simRunDf$deathRate <-  .01 #runif(nSims, 0, 1) #c(.1,.2)
  simRunDf$tweetRate <-  .25 #runif(nSims, 0, 2) #c(.2,.8)
  simRunDf$daysToDeath <-  20 #runif(nSims, 1, 20) #c(1,10)
  simRunDf$sdBeta <-  -1 #rep(.00001, nSims) #c(.0001,.2)
  simRunDf$nPatientZero <-  10 #c(10,100)
  simRunDf$nDailyContactsInf <- 10
  betaInfRatePerDayList <- list(rep(simRunDf$betaMean, nDays))
  simRunDf$betaInfRatePerDay <- rep(betaInfRatePerDayList, nSims)
  
  # Run config params, for each runDf set, alter config and add
  simRunDf$description <- 'sim, block ode with twitter'
  simRunDf$odeSolver <- 'block'
  simRunDf$runTwitter <- 1
  simRunDf$modelToRun <- 'twitter_sird'
  
  runDf = simRunDf
  i = 1
  reduction = 10000
  simDf = sirtd_vary_beta(seed = runDf[i,]$seed,
                        n_pop = round(runDf[i,]$nPop/reduction), 
                        n_days = runDf[i,]$nDays, 
                        print = FALSE,
                        beta_daily_inf_rates = 
                          unlist(runDf[i,]$betaInfRatePerDay),
                        gamma_res_per_day_rate = runDf[i,]$gamma,
                        tweet_rate_infected = runDf[i,]$tweetRate,
                        mean_days_to_death_from_t = runDf[i,]$daysToDeath,
                        n_patient_zero = runDf[i,]$nPatientZero,
                        death_prob = runDf[i,]$deathRate)
  expect_equal(sum(simDf$tweets), 928)
})

test_that("NA test", {
  seed = 4614
  nPop <- 10000
  nWeeks <- 7
  nDays <- 70
  gammaResolvedPerDayRateSim <- 0.38
  tweetRateInfected <- 1.99
  deathRateSim <- 0.54
  meanDaysToDeathFromT <- 14.79
  sdOfBetas <- 0.00
  betaInfectionRateMeanSim <- 0.18
  nPatientZero <- 10
  nDailyContacts <- 10  
  betaForWeek = abs(rnorm(nWeeks, betaInfectionRateMeanSim, sdOfBetas))
  betaDailyInfectionRatesSim = rep(betaForWeek, times = rep(7,nWeeks))
  
  simDf = sirtd_vary_beta(seed = seed,
                        n_pop = nPop,
                        n_days = nDays,
                        print = FALSE,
                        beta_daily_inf_rates = 
                          betaDailyInfectionRatesSim,
                        gamma_res_per_day_rate = gammaResolvedPerDayRateSim,
                        tweet_rate_infected = tweetRateInfected,
                        mean_days_to_death_from_t = meanDaysToDeathFromT, 
                        n_patient_zero = nPatientZero,
                        death_prob = deathRateSim)
  
  expect_equal(nrow(simDf), 70)
})


test_that("SIRTD test1",{
  seed = 43614
  nPop = 10000
  nWeeks = 10
  nDays = nWeeks * 7
  gammaResolvedPerDayRateSim = 1/7
  tweetRateInfected = .5
  deathRateSim = 0.1
  meanDaysToDeathFromT = 1
  sdOfBetas = .00001
  betaInfectionRateMeanSim = .3
  betaForWeek = abs(rnorm(nWeeks, betaInfectionRateMeanSim, sdOfBetas))
  betaDailyInfectionRatesSim = rep(betaForWeek, times = rep(7,nWeeks))
  nPatientZero = 10
  nDailyContacts = 10
  
  
  simDf = sirtd_vary_beta(seed = seed,
                        n_pop = nPop,
                        n_days = nDays,
                        print = FALSE,
                        beta_daily_inf_rates = 
                          betaDailyInfectionRatesSim,
                        gamma_res_per_day_rate = gammaResolvedPerDayRateSim,
                        tweet_rate_infected = tweetRateInfected,
                        mean_days_to_death_from_t = meanDaysToDeathFromT, 
                        n_patient_zero = nPatientZero,
                        death_prob = deathRateSim)
  
  
  expect_equal(nrow(simDf), 70)
  expect_equal(simDf[50,]$tweets, 857)
  expect_equal(simDf[50,]$d, 299)
})



