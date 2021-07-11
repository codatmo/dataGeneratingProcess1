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
    use_twitter = 1
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
fit_sim$summary(c("beta", "omega", "dI", "dT", "proportion_twitter"))

# # A tibble: 5 x 10
#   variable               mean   median       sd      mad       q5      q95  rhat ess_bulk ess_tail
#   <chr>                 <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl> <dbl>    <dbl>    <dbl>
# 1 beta                0.193    0.193   0.00208  0.00212   0.189    0.196    1.00    1687.    2336.
# 2 omega               0.00647  0.00646 0.000603 0.000593  0.00550  0.00749  1.00    2190.    2397.
# 3 dI                 11.6     11.6     0.397    0.404    11.0     12.2      1.00    1493.    1818.
# 4 dT                 16.8     16.8     1.90     1.87     13.8     20.0      1.00    2043.    2253.
# 5 proportion_twitter  0.418    0.417   0.0199   0.0201    0.386    0.451    1.00    1709.    2179.
