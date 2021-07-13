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
    beta_regularization = 0.2,
    prior_omega_mean = 0.4,
    prior_omega_std = 0.5,
    prior_dI_mean = 10,
    prior_dI_std = 2,
    prior_dT_mean = 16,
    prior_dT_std = 2,
    prior_twitter_lambda = 1.5
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
fit_sim$summary(c("omega", "dI", "dT", "twitter_rate"))

# # A tibble: 4 x 10
#   variable        mean  median      sd     mad       q5     q95  rhat ess_bulk ess_tail
#   <chr>          <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl> <dbl>    <dbl>    <dbl>
# 1 omega         0.0115  0.0113 0.00176 0.00157  0.00904  0.0147  1.00    2542.    2391.
# 2 dI           10.3    10.3    0.200   0.187   10.0     10.7     1.00     922.    1147.
# 3 dT           16.9    16.9    1.86    1.81    13.9     20.0     1.00    5570.    3032.
# 4 twitter_rate  0.546   0.535  0.0431  0.0333   0.498    0.635   1.00    1459.    1060.

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
#     week    mean    sd ground_truth
#    <dbl>   <dbl> <dbl>        <dbl>
#  1     0 0.135      NA       0.188
#  2     1 0.171      NA       0.254
#  3     2 0.346      NA       0.612
#  4     3 0.328      NA       0.314
#  5     4 0.315      NA       0.326
#  6     5 0.718      NA       0.643
#  7     6 0.922      NA       0.392
#  8     7 0.0137     NA       0.0470
#  9     8 0.00847    NA       0.163
# 10     9 0.00852    NA       0.211
