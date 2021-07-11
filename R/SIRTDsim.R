sirtd_vary_beta <- function(seed,
                            n_pop,
                            n_days,
                            print,
                            beta_daily_inf_rates,
                            gamma_res_per_day_rate,
                            death_prob,
                            tweet_rate_infected,
                            n_patient_zero,
                            mean_days_to_death_from_t) {
  set.seed(seed)
  day_state = c(rep('i',n_patient_zero), rep('s', n_pop - n_patient_zero))
  tweets = rep(0, n_pop)
  col_names = c('day', 's', 'i', 'r', 't', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(col_names)))
  colnames(df) = col_names
  next_day_state = day_state
  for (day in 1:n_days) {
    df[day,] = rep(NA,length(col_names)) #setup data, will cause errors if I miss something
    df[day,]$s=  length(subset(day_state, day_state == 's'))
    df[day,]$i = length(subset(day_state, day_state == 'i'))
    df[day,]$r = length(subset(day_state, day_state == 'r'))
    df[day,]$t = length(subset(day_state, day_state == 't'))
    df[day,]$d = length(subset(day_state, day_state == 'd'))
    df[day,]$day = day
    df[day,]$tweets = sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, recovered=%d, terminal=%d,
           dead=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$t,
          df[day,]$d, df[day,]$tweets,
          beta_daily_inf_rates[day]/gamma_res_per_day_rate))
    }
    tweets = rep(0, n_pop) #start fresh every day, certainly wrong.
    for (per in 1:n_pop) {
      #end infectious period
      if (day_state[per] == 'i') {
        tweets[per] = rpois(1, tweet_rate_infected)
        if (rbinom(n = 1, size = 1, prob = gamma_res_per_day_rate) == 1) {
          if (rbinom(n = 1, size = 1, prob = death_prob) == 1) {
            next_day_state[per] = 't'
          }
          else {
            next_day_state[per] = 'r'
          }
        }
      }
      if (day_state[per] == 't') {

        if (rbinom(n = 1, size =1, prob = 1/mean_days_to_death_from_t) == 1) {
          next_day_state[per] = 'd'
        }
      }
    }
    for (per in 1:n_pop) {
      if (day_state[per] == 'i') {
        for (other_per in sample(1:n_pop, rpois(n = 1, lambda = beta_daily_inf_rates[day]))) {
          if (day_state[other_per] == 's') {  #subtle but will reduce infs if there are not lots of s
            next_day_state[other_per] = 'i'
          }
        }
      }
    }
    day_state = next_day_state
  }
  return(df)
}

