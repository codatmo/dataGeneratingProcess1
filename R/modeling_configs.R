#' Returns SIRD configurations for model 'twitter_sird' that use tweets and 
#' not. ODE solver not specified
#' @param run_df The dataframe that is being copied and built up.
model_stan_baseline <- function(run_df) {

  model_no_tweets_df <- copy(run_df)
  model_no_tweets_df$apply_twitter_data <- 0
  model_no_tweets_df$description <- paste0(model_no_tweets_df$description,
                                           " no tweets")
  model_tweets_df <- copy(run_df)
  model_tweets_df$apply_twitter_data <- 1
  model_tweets_df$description <- paste0(model_tweets_df$description,
                                      " use tweets")

  combined_df = rbind(model_no_tweets_df, model_tweets_df)
  combined_df$ode_solver <- 'block'
  combined_df$model_to_run <- 'baseline'
  return(combined_df)
}


#' Returns SIRD configurations for model 'twitter_sird' that use tweets and 
#' not. ODE solver not specified
#' @param run_df The dataframe that is being copied and built up.
model_stan_UNINOVE_Brazil <- function(run_df) {
  
  model_df <- copy(run_df)
  model_df$apply_twitter_data <- 1
  model_df$model_to_run <- 'UNINOVE_Brazil'
  model_df$description <- paste0(model_df$description,
                                        " UNINOVE_Brazil")
  return(model_df)
}