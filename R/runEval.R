library(tidyverse)
library(cmdstanr)
rm(list = ls())
loadBrazil = FALSE
loadSim = TRUE
source(here::here("R","util.R"))
source(here::here("R","SIRTDsim.R"))
seed = 93435

nPop <- 210147125
nDays <- 311
nSims <- 1
nConfigs <- 1

# Set up our runs, shared across sim and real data runs
runDf <- data.frame(runId = NA)
runDf$description <- NA
runDf$tweetsEvalInInterval <- NA
runDf$casesEvalInInterval <- NA
runDf$odeSolver <- NA
runDf$runTwitter <- NA
runDf$nPop <- nPop
runDf$nDays <- nDays
runDf$seed <- 23454
runDf$modelToRun <- NA
runDf$compute_likelihood <- 1

# Simulator params
# replicate rows for nSims
nSims <- 1
simRunDf = runDf
simRunDf$betaMean <-  .3 #runif(nSims, 0, 1) #c(.2,.4)
simRunDf$gamma <- .29 #runif(nSims, 0, 1) #c(1/2,1/10)
simRunDf$deathRate <-  .01 #runif(nSims, 0, 1) #c(.1,.2)
simRunDf$tweetRate <-  .25 #runif(nSims, 0, 2) #c(.2,.8)
simRunDf$daysToDeath <-  20 #runif(nSims, 1, 20) #c(1,10)
simRunDf$sdBeta <-  -1 #rep(.00001, nSims) #c(.0001,.2)
simRunDf$nPatientZero <-  10 #c(10,100)
simRunDf$nDailyContactsInf <- 10
betaInfRatePerDayList <- list(rep(simRunDf$betaMean, nDays))
simRunDf$betaInfRatePerDay <- rep(betaInfRatePerDayList, nSims)

# Run config params, for each runDf set, alter config and add
simRunDf$description <- 'sim, block ode with twitter'
simRunDf$odeSolver <- 'block'
simRunDf$runTwitter <- 1
simRunDf$modelToRun <- 'twitter_sird'

runDf = simRunDf

dataDf = NA

if (loadBrazil) {
  tweets = read.csv(here::here("data","tweet_count.csv"))
  colnames(tweets) = c('dateT','predicted')
  brazilDeaths = readRDS(here::here("data","brazil_nation.rds"))

  dataDf = data.frame(date = brazilDeaths[1:354,]$date, 
                    d = brazilDeaths[1:354,]$last_available_deaths)
  tweetsPadded = rbind(data.frame(dateT = rep(NA,106), 
                                predicted = rep(0,106)), 
                     tweets)
  dataDf = cbind(dataDf,tweetsPadded)
  colnames(dataDf) = c('date','d','dateT','tweets')
  dataDf = dataDf %>% mutate(perc_d = d/nPop) %>% 
    mutate(perc_t = tweets/nPop)
  dataDf = dataDf[1:nDays,] #doing 2020 only
}

reduction = 10000
for (i in 1:nrow(runDf)) {
  if (loadSim) {
    simDf = SirtdVaryBeta(runName = runDf[i,]$runId, 
                          seed = runDf[i,]$seed,
                          nPop = runDf[i,]$nPop/reduction, 
                          nDays = runDf[i,]$nDays, 
                          print = TRUE,
                          betaDailyInfectionRates = unlist(runDf[i,]$betaInfRatePerDay),
                          gammaResolvedPerDayRate = runDf[i,]$gamma,
                          tweetRateInfected = runDf[i,]$tweetRate,
                          meanDaysToDeathFromT = runDf[i,]$daysToDeath,
                          nPatientZero = runDf[i,]$nPatientZero,
                          nDailyContacts = runDf[i,]$nDailyContactsInf,
                          deathRate = runDf[i,]$deathRate)
       dataDf = data.frame(d = simDf$d, 
                            tweets = simDf$tweets)
      }
  dataDf$day = 1:nrow(dataDf)
  fit <- NA
  if (runDf[i,]$modelToRun == 'twitter_sird') {
      stan_data <- 
          list(n_days = runDf[i,]$nDays,
               sDay1 = runDf$nPop - 1,
               iDay1 = 1,
               rDay1 = 0,
               dDay1 = 0,            
               NPop = runDf$nPop,
               tweets = dataDf$tweets,
               deaths = dataDf$d,
               compute_likelihood = runDf[i,]$compute_likelihood,
               run_twitter = runDf[i,]$runTwitter,
               run_block_ODE = ifelse(runDf[i,]$odeSolver == 'block', 1, 0),
               run_rk45_ODE = ifelse(runDf[i,]$odeSolver == 'rk45', 1, 0))

        model <- cmdstan_model(here::here("stan", "tweet_sird.stan"))

        fit <- model$sample(data=stan_data,
                     parallel_chains = 4,
                     iter_warmup = 1000,
                     iter_sampling = 1000,
                     chains = 4,
                     seed = 4857)
      }

      else if (runDf[i,]$modelToRun == 'tweet_sird_negbin_optimized') {
          stan_data_2 <- list(n_days = nrow(df),
                        y0 = c(runDf[i,]$nPop - runDf[i,]$nPatientZero, runDf[i,]$nPatientZero, 0, 0, 0), # one more zero here
                        t0 = 0,
                        ts = 1:runDf[i,]$nDays,
                        compute_likelihood = runDf[i,]$compute_likelihood,
                        use_twitter = runDf[i,]$run_twitter,
                        death_count = df$d,
                        symptomaticTweets = df$tweets,
                        trapezoidal_solver = 0)

          model2 <- cmdstan_model(here::here("stan", "tweet_sird_negbin_optimized.stan"))

          fit <- model2$sample(data=stan_data_2,
                          parallel_chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          chains = 4,
                          seed = 4857)
      }
    else {
      print(sprintf("no model selected, got:",runDf[i,]$modelToRun));
  }

  if (loadSim && TRUE) {
    resultsOdeDf = fit$summary(variables = c('ode_states'))
    compartmentNames = c('s','i','r','d');
    tweetSourceName = compartmentNames[2]
    odeDf = data.frame(matrix(resultsOdeDf$median, nrow = nrow(simDf),
                              ncol = length(compartmentNames)))
    colnames(odeDf) = compartmentNames
    odeDf$day = 1:nrow(simDf)
    
    odeLongDf = gather(data = odeDf, key = "compartmentODE", value = "mean",
                       all_of(compartmentNames))
    
    simLongDf = gather(data = simDf, key = "compartmentSim", value = "count",
                       all_of(c('tweets', compartmentNames)))
    
    plot = ggplot(data = NULL, aes(x = day, y = nPop * mean)) +
      geom_line(data = odeLongDf, aes(color = compartmentODE)) +
      geom_line(data = simLongDf, linetype = 'dotted', aes(y = count * reduction, 
                                       color = compartmentSim)) +
      labs(y = "median with sim data in dots",
           caption = paste0("predicted shaded, lines with dots for simulated truth, compartment is = ",
                            compartmentName), " tweets from ", tweetSourceName) +
      theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
    
    print(plot)
  }
  
  if (loadSim && FALSE) {
    recovPars = 
      fit$summary(variables = c('gamma', 'beta', 'deathRate', 'lambda_twitter'), 
                  mean, sd)
    cat(sprintf("\nBeta sim = %.4f vs recovered %.4f, sd=%.4f",
                runDf[i,]$betaMean, 
                recovPars[recovPars$variable == 'beta',]$mean,
                recovPars[recovPars$variable == 'beta',]$sd))
    
    cat(sprintf("\nGamma sim = %.4f vs recovered %.4f, sd=%.4f",
        runDf[i,]$gamma, 
        recovPars[recovPars$variable == 'gamma',]$mean,
        recovPars[recovPars$variable == 'gamma',]$sd))

    cat(sprintf("\nDeaths sim = %.4f vs recovered %.4f, sd=%.4f",
                runDf[i,]$deathRate, 
                recovPars[recovPars$variable == 'deaths',]$mean,
                recovPars[recovPars$variable == 'deaths',]$sd))
    cat(sprintf("\nLambda Twitter sim = %.4f vs recovered %.4f, sd=%.4f",
                runDf[i,]$tweetRate, 
                recovPars[recovPars$variable == 'lambda_twitter',]$mean,
                recovPars[recovPars$variable == 'lambda_twitter',]$sd))
  }
  
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
  
  runDf[i,]$casesInInterval = countPredInInterval(truth = dataDf$perc_d,
                                                  fitPredDf = predCasesDf,
                                                  maxQuantileLabel = maxQuantileLabel,
                                                  minQuantileLabel = minQuantileLabel)
  
  runDf[i,]$tweetsInInterval = countPredInInterval(truth = dataDf$perc_t,
                                                   fitPredDf = predTweetsDf,
                                                   maxQuantileLabel = maxQuantileLabel,
                                                   minQuantileLabel = minQuantileLabel)
  
      predCasesLine = 'predicted cases'
      predCasesRibbon = paste0('predicted cases ',minQuantileLabel, '/',
                               maxQuantileLabel)
      print(sprintf("Run %d predicted %d cases and %d tweets within interval",
            i, runDf[i,]$casesInInterval, runDf[i,]$tweetsInInterval))

      if (FALSE) {
        
      
        
        
        plot <- ggplot(data = NULL, aes(x = day, y = mean)) +
        # geom_ribbon(data = predCasesDf, aes(ymin = .data[[minQuantileLabel]],
        #                                     ymax = .data[[maxQuantileLabel]],
        #                                     fill = predCasesRibbon),
        #             alpha = 0.3, fill = "red") +
        geom_line(data = predCasesDf,
                   aes(color = 'pred cases death')) +
        geom_point(data = dataDf, 
                   aes(y = perc_d, color = 'sim cases'), size = .5, color = "black")
          
        #   
        #           #cases twitter informed
        # geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]],
        #                                       ymax = .data[[maxQuantileLabel]],
        #                                       fill = predCasesRibbon),
        #             alpha = 0.3, fill = "blue") +
        # 
        # 
        # #twitter
        # geom_ribbon(data = predTweetsDf,  aes(ymin = .data[[minQuantileLabel]],
        #                                       ymax = .data[[maxQuantileLabel]],
        #                                       fill = 'pred tweets'),
        #             alpha = 0.3, fill = "yellow") +
        # geom_line(data = predTweetsDf, aes(color = 'pred tweets')) +
        # geom_point(data = dataDf, aes(y = tweets, color = 'sim tweets'), size = .5, color = 'black') +
        # 
        # #  scale_fill_manual(labels = c('infected', 'tweets'))
        # labs(y = "mean with 80% shaded interval and sim data in dots",
        #      caption = summary) +
        # theme(plot.caption=element_text(size=12, hjust=0, margin=margin(15,0,0,0)))
        print(plot)
      }
      
  saveRDS(runDf,here::here("R",sprintf("%d_of_%devalBrazil622.rds",i,nrow(runDf))))
}
