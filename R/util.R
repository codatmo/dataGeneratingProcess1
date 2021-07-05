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
