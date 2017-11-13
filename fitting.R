require(CircStats) # for von Mises distribution
require(boot) # for logit
require(MASS)
require(reshape2)
require(stats)
require(moveHMM)


setwd("C:/Users/Patrick/Desktop/Whale")
whaledata <- read.csv("whale_data_cleaned.csv")

summary(whaledivestats$step)
summary(whaledivestats$angle)

mu0 <- c(0.015,0.08) # step mean (two parameters: one for each state)
sigma0 <- c(0.01,2) # step SD
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi/2,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)
## call to fitting function, but until now it doens not work
m <- fitHMM(data=d,nbStates=2,stepPar0=stepPar0,
            anglePar0=anglePar0)
