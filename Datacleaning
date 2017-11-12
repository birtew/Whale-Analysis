#function which computes the time in seconds in our case

overallsec <- function(datetime){
    data <- colsplit(datetime, " ", c("date", "time")) #split datetime in date und time
    date <- data$date #save date and time independently
    time <- data$time
    data <- colsplit(date, "-",c("year", "month", "day")) #split date
    year <- data$year
    month <- data$month
    day <- data2$day
    data <- colsplit(time, ":",c("hour", "min", "sec")) #split time
    hour <- data$hour
    min <- data$min
    sec <- data$sec
    seconds <- sec + min * 60 + hour * 60 *60 +day * 24 *60*60 #compute the date in seconds, since month and year are the same for all datapoints we can forget about them
    return(seconds)
}
