library(testthat)

source(here::here('R','sim_configs.R'))
source(here::here('R','util.R'))


test_that("sim_draw_params test", {
  run_df <- setup_run_df(seed = 123, n_pop = 100000, n_days = 40)
  brazil_df <- sim_brazil_1(run_df)
  draws_run_df <- sim_draw_params(n_sim = 2, source_df = brazil_df)
    expect_equal(draws_run_df[1,]$gamma != draws_run_df[2,]$gamma, TRUE)
})

test_that("Brazil 1 test", {
   run_df <- setup_run_df(seed = 123, n_pop = 100000, n_days = 40)
   brazil_df <- sim_brazil_1(run_df)
   expect_equal(brazil_df$seed, 123) 
})