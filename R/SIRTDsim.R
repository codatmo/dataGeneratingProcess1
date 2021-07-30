sirtd_vary_beta2 <- function(seed,
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
  s <- (n_patient_zero + 1):n_pop
  i <- 1:n_patient_zero
  r <- c()
  t <- c()
  d <- c()
  tweets = 0
  col_names = c('day', 's', 'i', 'r', 't', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(col_names)))
  colnames(df) = col_names
  for (day in 1:n_days) {
    df[day,] = rep(NA,length(col_names)) #setup data, will cause errors if I miss something
    df[day,]$s=  length(s)
    df[day,]$i = length(i)
     df[day,]$r = length(r)
    df[day,]$t = length(t)
    df[day,]$d = length(d)
    df[day,]$day = day
    df[day,]$tweets = tweets
    i_next <- c()
    t_next <- c()
    if (print) {
      cat(
        sprintf(
          "T Day=%d, susceptible=%d, infected=%d, recovered=%d, terminal=%d,
           dead=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$t,
          df[day,]$d, df[day,]$tweets,
          beta_daily_inf_rates[day]/gamma_res_per_day_rate))
    }
    tweets = 0 #start fresh every day, certainly wrong.
    for (per in i) {
      #end infectious period
      tweets = tweets + rpois(1, tweet_rate_infected)
      if (rbinom(n = 1, size = 1, prob = gamma_res_per_day_rate) == 1) {
        if (rbinom(n = 1, size = 1, prob = death_prob) == 1) {
          t_next <- c(t_next, per)
        }
        else {
          r <- c(r, per)
        }
      }
      else {
        i_next <- c(i_next, per)
      }
    }
    for(per in t) {
      if (rbinom(n = 1, size = 1, prob = 1/mean_days_to_death_from_t) == 1) {
        d <- c(d, per)
      }
      else {
        t_next <- c(t_next, per)
      }
    }
    s_delta <- 0
    for (per in i) {
#      print("I")
#      print(i)
      for (other_per in 1:rpois(n = 1, lambda = beta_daily_inf_rates[day])) {
        if (length(s)- s_delta >= 1) {
 #         print("S_next")
#          print(s_next)
          # i_next <- c(i_next, s[1 + s_delta]) # For now everyone gets a unique int id
          s_delta <- s_delta + 1
        }
      }
    }
    if (s_delta > 0) {
      i_next <- c(i, s[1:s_delta])
      s <- s[(s_delta + 1):length(s)]
    }
    i <- i_next
    t <- t_next
  }
  return(df)
}


sirtd_vary_beta_orig <- function(seed,
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

sirtd_vary_beta_1 <- function(seed,
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
  i <- n_patient_zero
  r <- 0
  t <- 0
  d <- 0
  s <- n_pop - n_patient_zero 
  tweets = 0
  col_names = c('day', 's', 'i', 'r', 't', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(col_names)))
  colnames(df) = col_names
  for (day in 1:n_days) {
    df[day,] = rep(NA,length(col_names)) 
    df[day,]$s=  s
    df[day,]$i = i
    df[day,]$r = r
    df[day,]$t = t
    df[day,]$d = d
    df[day,]$day = day
    df[day,]$tweets = tweets
    if (print) {
      cat(
        sprintf(
          "\nDay=%d, susceptible=%d, infected=%d, recovered=%d, terminal=%d,
           dead=%d, tweets=%d, R0=%.2f %d",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$t,
          df[day,]$d, df[day,]$tweets,
          beta_daily_inf_rates[day]/gamma_res_per_day_rate,
          s + i + r + t + d))
    }
    tweets = 0 #start fresh every day, certainly wrong.
    counter <- 0
    repeat { #move from t to d
      if (counter >= t) {
        break
      }
      if (rbinom(n = 1, size = 1, prob = 1/mean_days_to_death_from_t) == 1) {
        t <- t - 1
        d <- d + 1
      }
      counter <- counter + 1
    }
    counter <- 0
    i_today <- i
    repeat {
      if (counter >= i_today){
        break
      }
      #end infectious period
      tweets = tweets + rpois(1, tweet_rate_infected)
      #infect others
      for (other_per in 1:rpois(n = 1, lambda = beta_daily_inf_rates[day])) {
        if (s > 1) {
          s <- s - 1
          i <- i + 1
        }
      }
      #move from i to t,r 
      if (rbinom(n = 1, size = 1, prob = gamma_res_per_day_rate) == 1) {
        i <- i - 1
        if (rbinom(n = 1, size = 1, prob = death_prob) == 1) {
          t <- t + 1
        }
        else {
          r <- r + 1
        }
      }
      counter <- counter + 1
    }
  }
  return(df)
}

# sirtd_vary_beta_exact
#' Runs a simple deterministic SIRTD model with possible daily varied beta values. 
#' @param seed Seed for stochastic processes, ignored since this is a deterministic function
#' @param n_pop Population size
#' @param n_days Number of days to run simulation
#' @param print Boolean to control whether daily summary is printed
#' @param beta_daily_inf_rates Controls number of s that are infected by i to 
#' become i. Should be between 0 and 1
#' @param gamma_res_per_day_rate Controls number of i that transition to r or t
#' per day. Should be between 0 and 1
#' @param death_prob Probability transitioning from i to t
#' @param tweet_rate_infected If infected, what rate of tweets per day. Can be 
#' > 1
#' @param n_patient_zero Number of patients infected at day 0
#' @param mean_days_to_death_from_t Mean number of days in t before d
#' @return dataframe simulating sirtd+tweets state for number days 
sirtd_vary_beta_exact <- function(seed,
                            n_pop,
                            n_days,
                            print,
                            beta_daily_inf_rates,
                            gamma_res_per_day_rate,
                            death_prob,
                            tweet_rate_infected,
                            n_patient_zero,
                            mean_days_to_death_from_t) {
  old_seed <- .Random.seed
  on.exit({.Random.seed <<- old_seed})
  set.seed(seed)
  i <- n_patient_zero
  r <- 0
  t <- 0
  d <- 0
  s <- n_pop - n_patient_zero 
  col_names = c('day', 's', 'i', 'r', 't', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(col_names)))
  colnames(df) = col_names
  tweets = 0
  for (day in 1:n_days) {
    df[day,] = rep(NA,length(col_names)) 
    df[day,]$s=  round(s)
    df[day,]$i = round(i)
    df[day,]$r = round(r)
    df[day,]$t = round(t)
    df[day,]$d = round(d)
    df[day,]$day = day
    df[day,]$tweets = round(tweets)
    if (print) {
      cat(
        sprintf(
          "\nDay=%d, susceptible=%.2f, infected=%.2f, recovered=%.2f, terminal=%.2f,
           dead=%.2f, tweets=%.2f, R0=%.2f %.2f",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$t,
          df[day,]$d, df[day,]$tweets,
          beta_daily_inf_rates[day]/gamma_res_per_day_rate,
          s + i + r + t + d))
    }
    resolved <- gamma_res_per_day_rate * i  
    terminal <- death_prob * resolved
    recovered <- resolved - terminal
    dead <- 1/mean_days_to_death_from_t * t
    tweets = tweet_rate_infected * i
    infected <- min(i * beta_daily_inf_rates[day], s)
    s <- s - infected
    i <- i - resolved + infected
    r <- r + recovered
    t <- t - dead + terminal
    d <- d + dead
  }
  return(df)
}
# sirtd_vary_beta_exact

