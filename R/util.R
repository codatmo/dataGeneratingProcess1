#' Count the number of truth data points contained in the quantile interval of
#' model fit
#' @param truth vector of true or simulated values
#' @param fitPredDf dataframe with quantiles specified
#' 
countPredInInterval <- function(truth = truth, fitPredDf = fitPredDf,
                                maxQuantileLabel = maxQuantileLabel,
                                minQuantileLabel = minQuantileLabel) {
  actualsCovered = 0
  for (day in 1:length(truth)) {
    if (truth[day] <= fitPredDf[day,][[maxQuantileLabel]] &&
        truth[day] >= fitPredDf[day,][[minQuantileLabel]]) {
      actualsCovered = actualsCovered + 1
    }
  }
  return(actualsCovered)
} 

#' Take a fit, a run_df and the run index and compute the number of predicted
#' d and tweets are in the .2 to .8 central interval
#' @param fit cmdstanR fit object
#' @param run_df a dataframe setup for running models
#' @param j the current row being run
#' @param print boolean to control whether function prints, defaults to FALSE
#' 
countPredictionsInQuantile <- function(fit, run_df, j, print = FALSE) {
    minQuantile = .2
    maxQuantile = .8
    minQuantileLabel = '20%'
    maxQuantileLabel = '80%'
  
    predCasesDf = fit$summary(variables = c('pred_deaths'), mean,
                            ~quantile(.x, probs = c(minQuantile, maxQuantile),
                                      na.rm = TRUE))  #set to FALSE when Jose fixes his model
    predCasesDf$day = 1:nrow(predCasesDf)
  
    predTweetsDf = fit$summary(variables = c('pred_tweets'), mean,
                             ~quantile(.x, probs = c(minQuantile, maxQuantile),
                                       na.rm = TRUE))  #set to FALSE when Jose fixes his model
    predTweetsDf$day = 1:nrow(predTweetsDf)
  
    deaths_in_interval = countPredInInterval(truth = unlist(run_df[j,]$d),
                                                  fitPredDf = predCasesDf,
                                                  maxQuantileLabel = maxQuantileLabel,
                                                  minQuantileLabel = minQuantileLabel)
  
    tweets_in_interval = countPredInInterval(truth = unlist(run_df[j,]$tweets),
                                                   fitPredDf = predTweetsDf,
                                                   maxQuantileLabel = maxQuantileLabel,
                                                   minQuantileLabel = minQuantileLabel)
  
    if (print) {
      predCasesLine = 'predicted cases'
      predCasesRibbon = paste0('predicted cases ',minQuantileLabel, '/',
                           maxQuantileLabel)
      print(sprintf("Predicted %d cases and %d tweets within interval",
                  deaths_in_interval, tweets_in_interval))
    }
    return(c(deaths_in_interval, tweets_in_interval))
  }

#' Graph internal state counts for simulation and corresponding tweets for a run
#' of the runEval framework. Returns a ggplot geom_point element with x = days
#' y = count.
#' @param data_df one row of the run_df with simulation data added
#' @param hide_s Boolean to control whether to hide the s or susceptible counts
graph_sim_data <- function(data_df, hide_s) {
    sim_df = data.frame(day = 1:data_df$n_days, 
                        tweets = unlist(data_df$tweets), 
                        s = unlist(data_df$s),
                        i = unlist(data_df$i),
                        r = unlist(data_df$r),
                        t = unlist(data_df$t),
                        d = unlist(data_df$d))

    compartment_names <- c('s', 'i', 'r', 't', 'd')
    if (hide_s) {
      compartment_names <- compartment_names[-1]
    }
    sim_long_df = gather(data = sim_df, key = "compartment_sim", value = "count",
                         all_of(c('tweets', compartment_names)))
    return(geom_point(data = sim_long_df, aes(y = count, 
                                              color = compartment_sim),
                      size = .5))
}

#' Graph daily ODE means from SIRTD model. Returns a ggplot geom_line element
#' @param data_df A row from run_df
#' @param fit A fit object returned by cmdstanR with ode_states defined
#' @param hide_s Boolean to control whether to hide susceptible counts
graph_ODE <- function(data_df, fit, hide_s) {
  results_ODE_df = fit$summary(variables = c('ode_states'))
  compartment_names = c('s', 'i', 'r', 't', 'd')
  if (hide_s) {
    compartment_names <- compartment_names[-1]
  }
  ODE_df = data.frame(matrix(results_ODE_df$median, nrow = data_df$n_days,
                            ncol = length(compartment_names)))
  colnames(ODE_df) = compartment_names
  ODE_df$day = 1:data_df$n_days
  ODE_long_df = gather(data = ODE_df, key = "compartment_ODE", value = "count",
                     all_of(compartment_names))
  return(geom_line(data = ODE_long_df, aes(color = compartment_ODE)))
}

#' returns ggplot geom_line with optional ribbon around 20/80 central interval
#' of predictions. 
#' @param plot ggplot that this is adding to
#' @param prediction_label the label for the prediciton in fit object
#' @param fit cmdstanR fit object
#' @param show_ribbon control whether to show ribbon, default = TRUE
plot_predictions <- function(plot, prediction_label, fit, show_ribbon = TRUE) {
    minQuantile = .2
    maxQuantile = .8
    minQuantileLabel = '20%'
    maxQuantileLabel = '80%'
    pred_df = fit$summary(variables = c(prediction_label), mean,
                            ~quantile(.x, probs = c(minQuantile, maxQuantile),
                                      na.rm = FALSE))  #set to FALSE when Jose fixes his model
    pred_df$day = 1:nrow(pred_df)
    pred_df$count = pred_df$mean
    plot = plot + geom_line(data = pred_df)
    if (show_ribbon) {
      plot <- plot + geom_ribbon(data = pred_df, aes(ymin = .data[[minQuantileLabel]],
                                              ymax = .data[[maxQuantileLabel]],
                                              fill = "blue"),
                          alpha = 0.3, fill = "blue")
    }
    return(plot)
}

#' Returns string with calulation of parameter recovery
#' @param data_df row of run_df
#' @param fit cmdStanR fit object
param_recovery <- function(data_df, fit) {
  recov_pars = 
    fit$summary(variables = c('gamma', 'beta', 'deathRate', 'lambda_twitter'), 
                mean, sd)
    return(paste(sprintf("\nBeta sim = %.4f vs recovered %.4f, sd=%.4f",
            data_df$beta_mean, 
            recov_pars[recov_pars$variable == 'beta',]$mean,
            recov_pars[recov_pars$variable == 'beta',]$sd),
    sprintf("\nGamma sim = %.4f vs recovered %.4f, sd=%.4f",
            data_df$gamma, 
            recov_pars[recov_pars$variable == 'gamma',]$mean,
            recov_pars[recov_pars$variable == 'gamma',]$sd),
    sprintf("\nDeaths sim = %.4f vs recovered %.4f, sd=%.4f",
            data_df$death_prob, 
            recov_pars[recov_pars$variable == 'deathRate',]$mean,
            recov_pars[recov_pars$variable == 'deathRate',]$sd),
    sprintf("\nLambda Twitter sim = %.4f vs recovered %.4f, sd=%.4f",
            data_df$tweet_rate, 
            recov_pars[recov_pars$variable == 'lambda_twitter',]$mean,
            recov_pars[recov_pars$variable == 'lambda_twitter',]$sd),
    "\n"))
}

#' Generates the starting template for the dataframe that runs the experiments
#' @param seed The random seed for random processes in R
#' @param n_pop Population size, should be constant across runs
#' @param n_days Number of days to run, assumed to be constant across runs
#' @return dataframe appropriate for use in runEval.R
setup_run_df <- function(seed, n_pop, n_days) {
  # Set up our template, all runs need this info
  # Each row of dataframe is a separate run
  template_df <- data.frame(sim_run_id = c(NA))
  template_df$description <- NA
  template_df$seed <- seed

  # any simulation/model run
  template_df$n_pop <- n_pop
  template_df$n_days <- n_days

  # any model setup
  template_df$model_to_run <- 'none'
  template_df$compute_likelihood <- NA

  #setup data columns 
  template_df$s <- list(c())
  template_df$i <- list(c())
  template_df$r <- list(c())
  template_df$t <- list(c())
  template_df$d <- list(c())
  template_df$tweets <- list(c())

  # setup prediction columns
  template_df$d_in_interval <- NA_integer_
  template_df$tweets_in_interval <- NA_integer_
  set.seed(template_df$seed)
  return(template_df)
}
