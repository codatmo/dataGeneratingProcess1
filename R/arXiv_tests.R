library(cmdstanr)
library(tidyverse)
library(here)

set.seed(123)

model <- cmdstan_model(here("stan", "tweet_sird_negbin_optimized.stan"))

# SIRTD Sim function
SIRTD_sim <- function(runName = runName,
                          nPop = nPop,
                          nDays = nDays,
                          print = PRINT,
                          betaDailyInfectionRates = betaDailyInfectionRates,
                          gammaResolvedPerDayRate = gammaResolvedPerDayRate,
                          deathRate = deathRate,
                          tweetRateInfected = tweetRateInfected,
                          nPatientZero = nPatientZero,
                          meanDaysToDeathFromT = meanDaysToDeathFromT,
                          nDailyContacts = nDailyContacts,
                          seed = 123) {
    set.seed(seed)
    dayState <- c(rep("i", nPatientZero), rep("s", nPop - nPatientZero))
    tweets <- rep(0, nPop)
    colNames <- c("runName", "day", "s", "i", "r", "t", "d_today", "d", "tweets")
    df <- data.frame(matrix(nrow = 0, ncol = length(colNames)))
    colnames(df) <- colNames
    nextDayState <- dayState
    for (day in 1:nDays) {
        df[day, ] <- rep(NA, length(colNames)) # setup data, will cause errors if I miss something
        df[day, ]$runName <- runName
        df[day, ]$s <- length(subset(dayState, dayState == "s"))
        df[day, ]$i <- length(subset(dayState, dayState == "i"))
        df[day, ]$r <- length(subset(dayState, dayState == "r"))
        df[day, ]$t <- length(subset(dayState, dayState == "t"))
        df[day, ]$d_today <- length(subset(dayState, dayState == "d_today"))
        df[day, ]$d <- length(subset(dayState, dayState == "d"))
        df[day, ]$day <- day
        df[day, ]$tweets <- sum(tweets)
        if (print) {
            cat(
                sprintf(
                    "Day=%d, susceptible=%d, infected=%d, recovered=%d, terminal=%d,
          dead today=%d, dead=%d, tweets=%d, R0=%.2f\n",
                    df[day, ]$day, df[day, ]$s, df[day, ]$i, df[day, ]$r, df[day, ]$t, df[day, ]$d_today,
                    df[day, ]$d, df[day, ]$tweets, betaDailyInfectionRates[day] / gammaResolvedPerDayRate
                )
            )
        }
        tweets <- rep(0, nPop) # start fresh every day, certainly wrong.
        for (per in 1:nPop) {
            # end infectious period
            if (dayState[per] == "i") {
                tweets[per] <- rbinom(1, 1, tweetRateInfected)
                if (rbinom(n = 1, size = 1, prob = gammaResolvedPerDayRate) == 1) {
                    if (rbinom(n = 1, size = 1, prob = deathRate) == 1) {
                        nextDayState[per] <- "t"
                    }
                    else {
                        nextDayState[per] <- "r"
                    }
                }
            }
            if (dayState[per] == "t") {
                if (rbinom(n = 1, size = 1, prob = 1 / meanDaysToDeathFromT) == 1) {
                    nextDayState[per] <- "d"
                }
            }
        }
        for (per in 1:nPop) {
            if (dayState[per] == "i") {
                for (otherPer in sample(1:nPop, nDailyContacts)) {
                    if (dayState[otherPer] == "s" &&
                        rbinom(
                            n = 1, size = 1,
                            prob = min(betaDailyInfectionRates[day] / nDailyContacts, 1.0)
                        ) == 1) {
                        nextDayState[otherPer] <- "i"
                    }
                }
            }
        }
        dayState <- nextDayState
    }
    return(as_tibble(df))
}


##### Run 1 - Simulated Data #####
nPop <- 10000
nWeeks <- 10
nDays <- nWeeks * 7
gammaResolvedPerDayRateSim <- 1 / 7
tweetRateInfected <- .2
deathRateSim <- 0.1
meanDaysToDeathFromT <- 10

sdOfBetas <- 1e-6 # for time varying model of infectiousness
betaInfectionRateMeanSim <- .3
betaForWeek <- abs(rnorm(nWeeks, betaInfectionRateMeanSim, sdOfBetas))
betaDailyInfectionRatesSim <- rep(betaForWeek, times = rep(7, nWeeks))
nPatientZero <- 10
nDailyContacts <- 10

sim_df <- SIRTD_sim(
    runName = "test", nPop = nPop, nDays = nDays, print = TRUE,
    betaDailyInfectionRates = betaDailyInfectionRatesSim,
    gammaResolvedPerDayRate = gammaResolvedPerDayRateSim,
    tweetRateInfected = tweetRateInfected,
    meanDaysToDeathFromT = meanDaysToDeathFromT,
    nPatientZero = nPatientZero,
    nDailyContacts = nDailyContacts,
    deathRate = deathRateSim,
    seed = 123
)
# Getting new_deaths
# sim_df <- sim_df %>% mutate(d_today = coalesce(d - lag(d), 0))

stan_data_sim <- list(
    n_days = nrow(sim_df),
    y0 = c(first(sim_df$s), first(sim_df$i), first(sim_df$r), first(sim_df$t), first(sim_df$d)),
    t0 = 0,
    ts = seq_len(length(sim_df$day)),
    death_count = sim_df$d,
    symptomaticTweets = sim_df$tweets,
    compute_likelihood = 1,
    use_twitter = 1
)

fit_sim <- model$sample(
    data = stan_data_sim,
    output_dir = here("output", "simulated"),
    parallel_chains = 4,
    chains = 4, seed = 123
)

##### Run 2 Real Data Brazil #####
br_2020 <- read_csv(here("brazil_data", "COVID-2020.csv"))
tweets <- read_csv(here("brazil_data", "tweets_pred_arXiv_v1.csv"))


# Dropping dates that we don't have tweets
real_df <- tweets %>% inner_join(br_2020, "date")
population <- real_df %>%
    pull(estimated_population_2019) %>%
    max()
i_0 <- real_df %>%
    pull(new_confirmed) %>%
    first()
scaling_factor <- population / 1e5 # 100k
r_0 <- real_df$last_available_confirmed_per_100k_inhabitants[1] *  scaling_factor
d_0 <- real_df %>%
    pull(last_available_deaths) %>%
    first()
s_0 <- population - i_0 - r_0 - d_0
s_0 + i_0 + 0 + r_0 + d_0 == population

d_last <- real_df %>%
    pull(last_available_deaths) %>%
    last()
real_death_rate <- d_last / population

stan_data_brazil <- list(
    n_days = nrow(real_df),
    y0 = c(s_0, i_0, r_0, 0, d_0),
    t0 = 0,
    ts = seq_len(length(real_df$date)),
    death_count = real_df$last_available_deaths,
    symptomaticTweets = real_df$predicted,
    compute_likelihood = 1,
    use_twitter = 1
)

fit_brazil <- model$sample(
    data = stan_data_brazil,
    output_dir = here("output", "brazil"),
    parallel_chains = 4,
    chains = 4, seed = 123
)

#### Visualization #####
# Simulated Data
compartmentNames <- c("s", "i", "r", "t", "d")
sim_long_df <- gather(
    data = sim_df, key = "compartmentSim", value = "count",
    all_of(c("tweets", compartmentNames))
)

ggplot(data = NULL, aes(x = day, y = mean)) +
    geom_point(data = sim_long_df, aes(y = count, color = compartmentSim), size = .5) +
    labs(
        y = "sim data in dots",
        x = "sim days",
        caption = paste0("dots for simulated truth, compartment is = ") +
            theme(plot.caption = element_text(size = 12, hjust = 0, margin = margin(15, 0, 0, 0)))
    ) +
    scale_color_brewer(palette = "Set1")

state_S_sim <- fit_sim$summary("state_S")$mean
state_I_sim <- fit_sim$summary("state_I")$mean
state_R_sim <- fit_sim$summary("state_R")$mean
state_T_sim <- fit_sim$summary("state_T")$mean
state_D_sim <- fit_sim$summary("state_D")$mean
pred_deaths_sim <- fit_sim$summary("pred_deaths")$mean
pred_tweets_sim <- fit_sim$summary("pred_tweets")$mean

ode_df_sim <- tibble(
    day = seq_len(nrow(sim_df)),
    s = state_S_sim,
    i = state_I_sim,
    r = state_R_sim,
    t = state_T_sim,
    d = state_D_sim,
    tweets = pred_tweets_sim
)
ode_long_df_sim <- ode_df_sim %>% pivot_longer(-day, names_to = "compartmentODE", values_to = "mean")

plot_sim <- ggplot(data = NULL, aes(x = day, y = mean)) +
    geom_line(data = ode_long_df_sim, aes(color = compartmentODE)) +
    geom_point(data = sim_long_df, aes(y = count, color = compartmentSim), size = .5) +
    labs(
        y = "median with sim data in dots",
        caption = "predicted shaded, lines with dots for simulated truth",
        color = NULL
    ) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Set1")

ggsave(here("images", "fit_simulated.png"), plot_sim, dpi=300)

# Brazil Data
brazil_long_df <- tibble(
    day = seq_len(nrow(real_df)),
    d = real_df$last_available_deaths,
    tweets = real_df$predicted
) %>%
    pivot_longer(-day, names_to = "compartmentReal", values_to = "count")
state_S_br <- fit_brazil$summary("state_S")$mean
state_I_br <- fit_brazil$summary("state_I")$mean
state_R_br <- fit_brazil$summary("state_R")$mean
state_T_br <- fit_brazil$summary("state_T")$mean
state_D_br <- fit_brazil$summary("state_D")$mean
pred_deaths_br <- fit_brazil$summary("pred_deaths")$mean
pred_tweets_br <- fit_brazil$summary("pred_tweets")$mean

ode_df_br <- tibble(
    day = seq_len(nrow(real_df)),
    d = state_D_br,
    tweets = pred_tweets_br
)
ode_long_df_br <- ode_df_br %>% pivot_longer(-day, names_to = "compartmentODE", values_to = "mean")

plot_br <- ggplot(data = NULL, aes(x = day, y = mean)) +
    geom_line(data = ode_long_df_br, aes(color = compartmentODE)) +
    geom_point(data = brazil_long_df, aes(y = count, color = compartmentReal), size = .5) +
    labs(
        y = "median with sim data in dots",
        caption = "predicted shaded, lines with dots for real truth",
        color = NULL
    ) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Set1")

ggsave(here("images", "fit_brazil.png"), plot_br, dpi = 300)


#### Summary for Parameters ####

sim_vars <- c("beta", "omega", "proportion_twitter", "dI", "dT")
fit_sim_summary <- fit_sim$summary(sim_vars)
fit_sim_summary %>% write_csv(here("results", "sim_results.csv"))

brazil_vars <- c("beta", "omega", "proportion_twitter", "dI", "dT")
fit_brazil_summary <- fit_brazil$summary(brazil_vars)
fit_brazil_summary %>% write_csv(here("results", "brazil_results.csv"))
