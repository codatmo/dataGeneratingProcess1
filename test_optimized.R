library(cmdstanr)
library(tidyverse)
set.seed(43614)


# This is Breck's original function from June 10th
SIRDsim <- function(runName = runName,
                   nPop = nPop,
                   nDays = nDays,
                   print = PRINT,
                   betaInfRate = betaInfRate,
                   gammaResRate = gammaResRate,
                   deathRate = deathRate,
                   probTweet = probTweet,
                   infStartValue = infStartValue) {

  dayState <- c(rep("i", infStartValue), rep("s", nPop - infStartValue))
  tweets <- rep(0, nPop)
  contactPopulationSize <- 10
  colNames <- c("runName", "day", "s", "i", "r", "d_today", "d", "tweets")
  df <- data.frame(matrix(nrow = 0, ncol = length(colNames)))
  colnames(df) <- colNames
  nextDayState <- dayState
  for (day in 1:nDays) {
    df[day, ] <- rep(NA, length(colNames)) #setup data, will cause errors if I miss something
    df[day, ]$runName <- runName
    df[day, ]$s <- length(subset(dayState, dayState == "s"))
    df[day, ]$i <- length(subset(dayState, dayState == "i"))
    df[day, ]$r <- length(subset(dayState, dayState == "r"))
    df[day, ]$d_today <- length(subset(dayState, dayState == "d_today"))
    df[day, ]$d <- length(subset(dayState, dayState == "d"))
    df[day, ]$day <- day
    df[day, ]$tweets <- sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day = %d, susceptible = %d, infected = %d, resolved = %d, dead today = %d,
          dead=%d, tweets=%d, R0=%.2f\n",
          df[day, ]$day, df[day, ]$s, df[day, ]$i, df[day, ]$r, df[day, ]$d_today,
          df[day, ]$d, df[day, ]$tweets, betaInfRate / gammaResRate))
    }
    tweets <- rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == "i") {
        tweets[per] <- rbinom(1, 1, probTweet)
        if (rbinom(n = 1, size = 1, prob = gammaResRate) == 1) {
          if (rbinom(n = 1, size = 1, prob = deathRate) == 1) {
            nextDayState[per] <- "d" #not doing daily *1
          }
          else {
            nextDayState[per] <- "r"
          }
        }
      }
      if (dayState[per] == "d_today") { #irrelvant if not doing daily, see *1
        nextDayState[per] <- "d"
      }
    }
    for (per in 1:nPop) {
      if (dayState[per] == "i") {
        for (otherPer in sample(1:nPop, contactPopulationSize)) {
          if (dayState[otherPer] == "s" &&
              rbinom(n = 1, size = 1, prob = betaInfRate / contactPopulationSize) == 1) {
            nextDayState[otherPer] <- "i"
          }
        }
      }
    }
    dayState <- nextDayState
  }
  return(df)
}

nPop <- 10000
nWeeks <- 10
nDays <- nWeeks * 7
betaInfRateSim <- .3
gammaResRateSim <- 1 / 7
probTweetSim <- .5
deathRateSim <- .1
infStartValue <- 10

simDf <- SIRDsim(runName = "test", nPop = nPop, nDays = nDays, print = TRUE,
                betaInfRate = betaInfRateSim, gammaResRate = gammaResRateSim,
                probTweet = probTweetSim, deathRate = deathRateSim,
                infStartValue = infStartValue)

simDf %>%
  ggplot(aes(x = day)) +
  geom_line(aes(y = d, color = "death")) +
  geom_line(aes(y = tweets, color = "tweets")) +
  theme_minimal() +
  theme(legend.position = "bottom")

model <- cmdstan_model(here::here("stan", "tweet_sird_negbin_optimized.stan"))

# Code modified from
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

stan_data <- list(
  n_days = nrow(simDf),
  y0 = c(nPop - infStartValue, infStartValue, 0, 0),
  t0 = 0,
  ts = 1:nDays,
  N = nPop,
  death_count = simDf$d,
  symptomaticTweets =  simDf$tweets,
  compute_likelihood = 1,
  use_twitter = 1)

fit <- model$sample(data = stan_data, output_dir = "output",
                    parallel_chains = 4,
                    chains = 4, seed = 4857)

resultsBeta <- fit$summary(variables = "beta")
resultsGamma <- fit$summary(variables = "gamma")
resultsDeath <- fit$summary(variables = "death_rate")

minQuantile <- .2
maxQuantile <- .8
minQuantileLabel <- "20%"
maxQuantileLabel <- "80%"

predCasesDf <- fit$summary(variables = "pred_cases", mean,
                           ~quantile(.x, probs = c(minQuantile, maxQuantile)))
predCasesDf$day <- seq_len(nrow(predCasesDf))
predCasesTwitterDf <- fit$summary(variables = "pred_cases", mean,
                                  ~quantile(.x, probs = c(minQuantile, maxQuantile)))
predCasesTwitterDf$day <- seq_len(nrow(predCasesTwitterDf))
predTweetsDf <- fit$summary(variables = "pred_tweets", mean,
                            ~quantile(.x, probs = c(minQuantile, maxQuantile)))
predTweetsDf$day <- seq_len(nrow(predTweetsDf))

predCasesLine <- "predicted cases"
predCasesRibbon <- paste0("predicted cases ", minQuantileLabel, "/", maxQuantileLabel)

ggplot(data = NULL, aes(x = day, y = mean)) +
  #cases no twitter info
  geom_ribbon(data = predCasesDf, aes(ymin = .data[[minQuantileLabel]],
                                      ymax = .data[[maxQuantileLabel]],
                                      fill = predCasesRibbon),
              alpha = 0.3, fill = "red") +
  geom_line(data = predCasesDf,
            aes(color = "predicted cases")) +
  #cases twitter informed
  geom_ribbon(data = predCasesTwitterDf,  aes(ymin = .data[[minQuantileLabel]],
                                              ymax = .data[[maxQuantileLabel]],
                                              fill = predCasesRibbon),
              alpha = 0.3, fill = "blue") +
  geom_line(data = predCasesTwitterDf, aes(color = "predicted cases + twitter")) +

  geom_point(data = simDf, aes(y = d, color = "sim cases"), size = .5, color = "red") +

  #twitter
  geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]],
                                        ymax = .data[[maxQuantileLabel]],
                                        fill = "pred tweets"),
              alpha = 0.3, fill = "yellow") +
  geom_line(data = predTweetsDf, aes(color = "predicted tweets")) +
  geom_point(data = simDf, aes(y = tweets, color = "sim tweets"), size = .5, color = "black") +

  #  scale_fill_manual(labels = c('infected', 'tweets'))
  labs(y = "mean with 80% shaded interval and sim data in dots",
       caption = "predicted shaded, lines with dots for simulated truth") +
  theme(plot.caption = element_text(size = 12, hjust = 0, margin = margin(15, 0, 0, 0)))

resultsOdeDf <- fit$summary(variables = "ode_states")
compartmentNames <- c("s", "i", "r", "d")
odeDf <- data.frame(matrix(resultsOdeDf$median, nrow = nrow(simDf),
       ncol = length(compartmentNames)))
colnames(odeDf) <- compartmentNames
odeDf$day <- seq_len(nrow(simDf))

odeLongDf <- gather(data = odeDf, key = "compartmentODE", value = "mean",
                   all_of(compartmentNames))

simLongDf <- gather(data = simDf, key = "compartmentSim", value = "count",
                   all_of(c("tweets", compartmentNames)))

ggplot(data = NULL, aes(x = day, y = mean)) +
  geom_line(data = odeLongDf, aes(color = compartmentODE)) +
  geom_point(data = simLongDf, aes(y = count, color = compartmentSim), size = .5) +
  labs(y = "median with sim data in dots",
     caption = "predicted shaded, lines with dots for simulated truth") +
  theme(legend.position = "bottom")
