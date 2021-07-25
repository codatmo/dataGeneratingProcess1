#' Returns SIRD configurations for model 'twitter_sird' that use tweets and 
#' not. ODE solver not specified
#' @param run_df The dataframe that is being copied and built up.
model_stan_SIRDs <- function(run_df) {
  
  model_no_tweets_df <- copy(run_df)
  model_no_tweets_df$apply_twitter_data <- 0
  model_no_tweets_df$model_to_run <- 'twitter_sird'
  model_no_tweets_df$description <- paste0(model_no_tweets_df$description,
                                           "no tweets")
  model_tweets_df <- copy(run_df)
  model_tweets_df$apply_twitter_data <- 1
  model_tweets_df$model_to_run <- 'twitter_sird'
  model_tweets_df$description <- paste0(model_tweets_df$description,
                                      "use tweets")
  return(rbind(model_no_tweets_df, model_tweets_df))
}