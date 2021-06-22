SIRDsim = function(runName = runName,
                   nPop = nPop,
                   nDays = nDays,
                   print = PRINT,
                   betaInfRate = betaInfRate,
                   gammaResRate = gammaResRate,
                   deathRate = deathRate,
                   probTweet = probTweet,
                   infStartValue = infStartValue) {
  
  dayState = c(rep('i',infStartValue), rep('s', nPop - infStartValue))
  tweets = rep(0, nPop)
  contactPopulationSize = 10
  colNames = c('runName', 'day', 's', 'i', 'r', 'd_today', 'd', 'tweets')
  df = data.frame(matrix(nrow = 0, ncol=length(colNames)))
  colnames(df) = colNames
  nextDayState = dayState
  for (day in 1:nDays) {
    df[day,] = rep(NA,length(colNames)) #setup data, will cause errors if I miss something
    df[day,]$runName = runName
    df[day,]$s=  length(subset(dayState, dayState == 's'))
    df[day,]$i = length(subset(dayState, dayState == 'i'))
    df[day,]$r = length(subset(dayState, dayState == 'r'))
    df[day,]$d_today = length(subset(dayState, dayState == 'd_today'))
    df[day,]$d = length(subset(dayState, dayState == 'd'))
    df[day,]$day = day
    df[day,]$tweets = sum(tweets)
    if (print) {
      cat(
        sprintf(
          "Day=%d, susceptible=%d, infected=%d, resolved=%d, dead today=%d,
          dead=%d, tweets=%d, R0=%.2f\n",
          df[day,]$day, df[day,]$s, df[day,]$i, df[day,]$r, df[day,]$d_today,
          df[day,]$d, df[day,]$tweets, betaInfRate/gammaResRate))
    }
    tweets = rep(0, nPop) #start fresh every day, certainly wrong.
    for (per in 1:nPop) {
      #end infectious period
      if (dayState[per] == 'i') {
        tweets[per] = rbinom(1, 1, probTweet)
        if (rbinom(n = 1, size = 1, prob = gammaResRate) == 1) {
          if (rbinom(n = 1, size = 1, prob = deathRate) == 1) {
            nextDayState[per] = 'd' #not doing daily *1
          }
          else {
            nextDayState[per] = 'r'
          }
        }
      }
      if (dayState[per] == 'd_today') { #irrelvant if not doing daily, see *1
        nextDayState[per] = 'd'
      }
    }
    for (per in 1:nPop) {
      if (dayState[per] == 'i') {
        for (otherPer in sample(1:nPop, contactPopulationSize)) {
          if (dayState[otherPer] == 's' &&
              rbinom(n = 1, size = 1, prob = betaInfRate/contactPopulationSize) == 1) {
            nextDayState[otherPer] = 'i'
          }
        }
      }
    }
    dayState = nextDayState
  }
  return(df)
}