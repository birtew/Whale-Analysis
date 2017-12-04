install.packages("CircStats")
install.packages("boot")
install.packages("MASS")
install.packages("reshape2")
install.packages("stats")
install.packages("moveHMM")

require(CircStats) # for von Mises distribution
require(boot) # for logit
require(MASS)
require(reshape2)
require(stats)
require(moveHMM)

setwd("~/Uni/(M.Sc.) 3. Semester/Statistical Consulting/Minke whale project")
whalegps <- read.table("GPS_data_170921.txt", header=T, sep= "\t") # read the gps data
whaledivestats <- read.csv("dive_stats_170921.csv")


#function which computes the time in seconds in our case for the gps data. datetime must be in the format year-month-day hour:minute:second

overallsec <- function(datetime){
  data <- colsplit(datetime, " ", c("date", "time")) #split datetime in date und time
  date <- data$date #save date and time independently
  time <- data$time
  data <- colsplit(date, "-",c("year", "month", "day")) #split date
  year <- data$year
  month <- data$month
  day <- data$day
  data <- colsplit(time, ":",c("hour", "min", "sec")) #split time
  hour <- data$hour
  min <- data$min
  sec <- data$sec
  seconds <- sec + min * 60 + hour * 60 *60 +day * 24 *60*60 #compute the date in seconds, since month and year are the same for all datapoints we can forget about them
  return(seconds)
}

whalegps$overallsec <- overallsec(whalegps$datetime) # create the overallsec computation to the gps data

#function which computes the time in seconds in our case for the stat data. datetime must be in the format year-month-day hour:minute:second UTC

overallsecwithUTC <- function(datetime){
  data <- colsplit(datetime, " ", c("date", "time", "UTC")) #split datetime in date und time and UTC
  date <- data$date #save date and time independently
  time <- data$time
  data <- colsplit(date, "-",c("year", "month", "day")) #split date
  year <- data$year
  month <- data$month
  day <- data$day
  data <- colsplit(time, ":",c("hour", "min", "sec")) #split time
  hour <- data$hour
  min <- data$min
  sec <- data$sec
  seconds <- sec + min * 60 + hour * 60 *60 +day * 24 *60*60 #compute the date in seconds, since month and year are the same for all datapoints we can forget about them
  return(seconds)
}

whaledivestats$enddescsec <- overallsecwithUTC(whaledivestats$enddesc)
whaledivestats$begdescsec <- overallsecwithUTC(whaledivestats$begdesc)

#plot the time against lat and long

plot(whalegps$overallsec,whalegps$lat)
plot(whalegps$overallsec,whalegps$long)

#seems that linear interpolation is ok or?

interpolationvalue <- function(tvec,time,argument,k){ #tvec is vector of timepoints where we want to have the interpolation value, time is the set of times there we know the value we are interested in, argument is the vector of the variable we are interested in and k is the maximal distance in the timevector we want to allow for linear interpolation
  b <-rep(0, length(tvec))
  for(i in 1:length(tvec)){
    t <- tvec[i]
    if(t<time[1] || t > time[length(time)]){a <- NA}
    else{
      timevector <- t-time
      postime <- subset(timevector, timevector >= 0)
      t1 <- length(postime) #the index of the time point in time which is the biggest lower threshold
      t2<- t1+1 # now t1<= t <= t2
      if(time[t2]-time[t1] <= k){
        a <- argument[t1]+(t-time[t1])/(time[t2]-time[t1]) * (argument[t2]-argument[t1])
      }
      else { a <- NA}
    }
    b[i] <- a
  }
  return(b)
}

# calculation for the lat and long variables by linear interpolation in the dive set
whaledivestats$lat <- interpolationvalue(whaledivestats$begdescsec, whalegps$overallsec,whalegps$lat,10000000000)
whaledivestats$long <- interpolationvalue(whaledivestats$begdescsec, whalegps$overallsec,whalegps$long,10000000000)


# calculation of step and angle with moveHMM
d <- prepData(data.frame(lat=whaledivestats$lat,long=whaledivestats$long), type="LL", coordNames = c("lat","long"))
whaledivestats$step <- d$step
whaledivestats$angle <- d$angle

plot(d)


#
splitset <- function(k, argument, dataset){
  a <- which(argument >= k)
  zr <- vector("list", length = length(a) + 1)
  zr[[1]] <- dataset[1:(a[1]-1),]
  #zr <- data.frame()
  for (i in 2:length(a)) {
    zr[[i]] <- dataset[(a[i-1]+1):(a[i]-1),]
  }
  zr[[length(a)+1]] <- dataset[(a[length(a)]+1):length(argument),]
  return(zr)
}


k=20000
argument <- whaledivestats$postdive.dur
newdata <- splitset(k, argument, dataset)

#a[0]
#i=1
data <- whaledivestats
dataset <- whaledivestats



## Nullen in postdive.dur löschen
nozeros <- whaledata[(whaledata$postdive.dur > 0),]




### create csv
#write to new csv-file to prevent repeating the steps above everytime
write.csv(data, file = "whale_data_cleaned.csv")
