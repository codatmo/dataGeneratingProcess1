library(cmdstanr)
library(here)
library(dplyr)
library(tibble)
source(here("R", "SIRTDsim.R"))

model <- cmdstan_model(here("stan", "tweet_sirtd_negbin_matrix_exp_varying_beta.stan"))

set.seed(123)
weekly_betas <- rnorm(10, 0.3, 0.2)

df_sim_weekly_beta <- sirtd_vary_beta(
    seed = 93435,
    n_pop = 20000,
    n_days = 70,
    print = TRUE,
    beta_daily_inf_rates = rep(weekly_betas, each = 7),
    gamma_res_per_day_rate = 0.1,
    tweet_rate_infected = 0.5,
    mean_days_to_death_from_t = 16,
    n_patient_zero = 20,
    death_prob = 0.01
)

stan_data <- list(
    n_days = nrow(df_sim_weekly_beta),
    y0 = c(
        first(df_sim_weekly_beta$s),
        first(df_sim_weekly_beta$i),
        first(df_sim_weekly_beta$r),
        first(df_sim_weekly_beta$t),
        first(df_sim_weekly_beta$d)
    ),
    t0 = 0,
    ts = seq_len(length(df_sim_weekly_beta$day)),
    death_count = df_sim_weekly_beta$d,
    symptomaticTweets = df_sim_weekly_beta$tweets,
    compute_likelihood = 1,
    use_twitter = 1,
    beta_regularization = 10
)

fit_sim <- model$sample(
    data = stan_data,
    # output_dir = here("output", "simulated"),
    parallel_chains = 4,
    chains = 4,
    seed = 123
)


# dI is 1/gamma_res_per_day_rate
# omega is death_prob
fit_sim$summary(c("omega", "dI", "dT", "proportion_twitter"))

# # A tibble: 4 x 10
#   variable              mean  median      sd     mad       q5     q95  rhat ess_bulk ess_tail
#   <chr>                <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl> <dbl>    <dbl>    <dbl>
# 1 omega               0.0109  0.0108 0.00145 0.00124  0.00881  0.0133  1.01    1627.    1843.
# 2 dI                 10.2    10.2    0.199   0.186    9.80    10.5     1.01     386.     138.
# 3 dT                 16.8    16.9    1.81    1.86    13.9     19.7     1.00    1266.    1729.
# 4 proportion_twitter  0.525   0.522  0.0164  0.0134   0.504    0.551   1.01     688.     944.

fit_sim$summary("beta") %>%
    rownames_to_column("day") %>%
    mutate(
        day = as.numeric(day),
        week = (day - 1) %/% 7
    ) %>%
    group_by(week) %>%
    summarise(
        mean = mean(mean),
        sd = sd(mean)
    ) %>%
        add_column(ground_truth = weekly_betas)

# # A tibble: 10 x 4
#     week   mean    sd ground_truth
#    <dbl>  <dbl> <dbl>        <dbl>
#  1     0 0.134     NA       0.188
#  2     1 0.169     NA       0.254
#  3     2 0.343     NA       0.612
#  4     3 0.332     NA       0.314
#  5     4 0.323     NA       0.326
#  6     5 0.786     NA       0.643
#  7     6 0.724     NA       0.392
#  8     7 0.0988    NA       0.0470
#  9     8 0.0484    NA       0.163
# 10     9 0.0475    NA       0.211
