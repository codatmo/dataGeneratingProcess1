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
