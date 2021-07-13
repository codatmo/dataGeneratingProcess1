library(cmdstanr)
library(here)
library(dplyr)
source(here("R", "SIRTDsim.R"))

model <- cmdstan_model(here("stan", "tweet_sirtd_negbin_ODE.stan"))

df_sim_fixed_beta <- sirtd_vary_beta(
    seed = 93435,
    n_pop = 20000,
    n_days = 70,
    print = TRUE,
    beta_daily_inf_rates = rep(0.2, 70),
    gamma_res_per_day_rate = 0.1,
    tweet_rate_infected = 0.5,
    mean_days_to_death_from_t = 16,
    n_patient_zero = 20,
    death_prob = 0.01
)

stan_data <- list(
    n_days = nrow(df_sim_fixed_beta),
    y0 = c(first(df_sim_fixed_beta$s),
           first(df_sim_fixed_beta$i),
           first(df_sim_fixed_beta$r),
           first(df_sim_fixed_beta$t),
           first(df_sim_fixed_beta$d)),
    t0 = 0,
    ts = seq_len(length(df_sim_fixed_beta$day)),
    death_count = df_sim_fixed_beta$d,
    symptomaticTweets = df_sim_fixed_beta$tweets,
    compute_likelihood = 1,
    use_twitter = 1,
    prior_beta_mean = 2,
    prior_beta_std = 1,
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
    #output_dir = here("output", "simulated"),
    parallel_chains = 4,
    chains = 4,
    seed = 123
)


# dI is 1/gamma_res_per_day_rate
# omega is death_prob
fit_sim$summary(c("beta", "omega", "dI", "dT", "twitter_rate"))

# A tibble: 5 x 10
#   variable         mean   median       sd      mad       q5      q95  rhat ess_bulk ess_tail
#   <chr>           <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl> <dbl>    <dbl>    <dbl>
# 1 beta          0.193    0.193   0.00206  0.00199   0.189    0.196    1.00    1946.    2310.
# 2 omega         0.00648  0.00646 0.000608 0.000613  0.00550  0.00751  1.00    2195.    2200.
# 3 dI           11.5     11.5     0.383    0.375    10.9     12.2      1.00    1648.    1880.
# 4 dT           16.9     16.9     1.90     1.89     13.8     20.0      1.00    2331.    2794.
# 5 twitter_rate  0.420    0.420   0.0194   0.0196    0.388    0.452    1.00    1769.    2151.
