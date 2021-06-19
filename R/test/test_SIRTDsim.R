library(testthat)

source(here::here('R','SIRTDsim.R'))

test_that("NA test", {
  set.seed(4614)
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

  simDf = SirtdVaryBeta(runName = 'test', nPop = nPop, nDays = nDays, print = TRUE,
                      betaDailyInfectionRates = betaDailyInfectionRatesSim,
                      gammaResolvedPerDayRate = gammaResolvedPerDayRateSim,
                      tweetRateInfected = tweetRateInfected,
                      meanDaysToDeathFromT = meanDaysToDeathFromT,
                      nPatientZero = nPatientZero,
                      nDailyContacts = nDailyContacts,
                      deathRate = deathRateSim)
  
  expect_equal(nrow(simDf), 70)
})


test_that("SIRTD test1",{
  set.seed(43614)
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

  simDf = SirtdVaryBeta(runName = 'test', nPop = nPop, nDays = nDays, print = FALSE,
                      betaDailyInfectionRates = betaDailyInfectionRatesSim,
                      gammaResolvedPerDayRate = gammaResolvedPerDayRateSim,
                      tweetRateInfected = tweetRateInfected,
                      meanDaysToDeathFromT = meanDaysToDeathFromT,
                      nPatientZero = nPatientZero,
                      nDailyContacts = nDailyContacts,
                      deathRate = deathRateSim)
  
  expect_equal(nrow(simDf), 70)
  expect_equal(simDf[50,]$tweets, 854)
  expect_equal(simDf[50,]$d, 364)
})




    

