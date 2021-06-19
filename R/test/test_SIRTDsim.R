library(testthat)

source(here::here('R','SIRTDsim.R'))

test_that("SIRTD test",{
  set.seed(43614)
  nPop = 10000
  nWeeks = 10
  nDays = nWeeks * 7
  gammaResolvedPerDayRateSim = 1/7
  tweetRateInfected = .5
  deathRateSim = 0.1
  meanDaysToDeathFromT = 1

  sdOfBetas = .00001 # for time varying model of infectiousness
  betaInfectionRateMeanSim = .3
  betaForWeek = abs(rnorm(nWeeks, betaInfectionRateMeanSim, sdOfBetas))
  betaDailyInfectionRatesSim = rep(betaForWeek, times = rep(7,nWeeks))
  nPatientZero = 10
  nDailyContacts = 10

  simDf = SirtdVaryBeta(runName = 'test', nPop = nPop, nDays = nDays, print = TRUE,
                      betaDailyInfectionRates = betaDailyInfectionRatesSim,
                      gammaResolvedPerDayRate = gammaResolvedPerDayRateSim,
                      tweetRateInfected = tweetRateInfected,
                      meanDaysToDeathFromT = meanDaysToDeathFromT,
                      nPatientZero = nPatientZero,
                      nDailyContacts = nDailyContacts,
                      deathRate = deathRateSim)
  
  expect_equal(nrow(simDf), 70)
  expect_equal(simDf[50,]$tweets, 811)
  expect_equal(simDf[50,]$d, 368)
})

