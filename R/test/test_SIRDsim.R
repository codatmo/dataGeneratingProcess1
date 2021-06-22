library(testthat)

source(here::here('R','SIRDsim.R'))

test_that("SIRD test",{
  set.seed(43614)
  nPop = 10000
  nWeeks = 10
  nDays = nWeeks * 7
  betaInfRateSim = .3
  gammaResRateSim = 1/7
  probTweetSim = .5
  deathRateSim = .1
  infStartValue = 10
  
  simDf = SIRDsim(runName = 'test', nPop = nPop, nDays = nDays, print = FALSE,
                  betaInfRate = betaInfRateSim, gammaResRate = gammaResRateSim,
                  probTweet = probTweetSim, deathRate = deathRateSim,
                  infStartValue = infStartValue)
  
  expect_equal(nrow(simDf), 70)
  expect_equal(simDf[50,]$tweets, 862)
  expect_equal(simDf[50,]$d, 362)
    expect_equal(simDf[70,]$tweets, 309)
  expect_equal(simDf[70,]$d, 714)
  
  })

