SirtdVaryBeta <- function(runName = runName,
                         nPop = nPop,
                         nDays = nDays,
                         print = PRINT,
                         betaDailyInfectionRates = betaDailyInfectionRates,
                         gammaResolvedPerDayRate = gammaResolvedPerDayRate,
                         deathRate = deathRate,
                         tweetRateInfected = tweetRateInfected,
                         nPatientZero = nPatientZero,
                         meanDaysToDeathFromT = meanDaysToDeathFromT,
                         nDailyContacts = nDailyContacts) {

  dayState = c(rep('i',nPatientZero), rep('s', nPop - nPatientZero))
  tweets = rep(0, nPop)
  colNames = c('runName', 'day', 's', 'i', 'r', 't', 'd_today', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(colNames)))
  colnames(df) = colNames
  nextDayState = dayState
  for (day in 1:nDays) {
    df[day,] = rep(NA,length(colNames)) #setup data, will cause errors if I miss something
    df[day,]$runName = runName
    df[day,]$s=  length(subset(dayState, dayState == 's'))
    df[day,]$i = length(subset(dayState, dayState == 'i'))
    df[day,]$r = length(subset(dayState, dayState == 'r'))
    df[day,]$t = length(subset(dayState, dayState == 't'))
    df[day,]$d_today = length(subset(dayState, dayState == 'd_today'))
    df[day,]$d = length(subset(dayState, dayState == 'd'))
    df[day,]$day = day
    df[day,]$tweets = sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, recovered=%d, terminal=%d,
          dead today=%d, dead=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$t, df[day,]$d_today,
          df[day,]$d, df[day,]$tweets, betaDailyInfectionRates[day]/gammaResolvedPerDayRate))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'i') {
        tweets[per] = as.integer(rbernoulli(1, tweetRateInfected))
        if (rbinom(n = 1, size = 1, prob = gammaResolvedPerDayRate) == 1) {
          if (rbinom(n = 1, size = 1, prob = deathRate) == 1) {
            nextDayState[per] = 't'
          }
          else {
            nextDayState[per] = 'r'
          }
        }
      }
      if (dayState[per] == 't') {
        if (rbinom(n = 1, size =1, prob = 1/meanDaysToDeathFromT) == 1) {
          nextDayState[per] = 'd'
        }
      }
    }
    for (per in 1:nPop) {
      if (dayState[per] == 'i') {
        for (otherPer in sample(1:nPop, nDailyContacts)) {
          if (dayState[otherPer] == 's' &&
              rbinom(n = 1, size = 1,
                     prob = min(betaDailyInfectionRates[day],1.0)) == 1) {
            nextDayState[otherPer] = 'i'
          }
        }
      }
    }
    dayState = nextDayState
  }
  return(df)
}

