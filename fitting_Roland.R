require(CircStats) # for von Mises distribution
require(boot) # for logit
require(MASS)
require(reshape2)
require(stats)
require(moveHMM)


#setwd("C:/Users/Patrick/Desktop/Whale")
setwd("~/Uni/(M.Sc.) 3. Semester/Statistical Consulting/Minke whale project")
whaledata <- read.csv("whale_data_cleaned.csv")

obs <- whaledata[(7:3583), c(5, 7, 8, 13, 14)]




## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
# mu & kappa: von Mises distr.
pn2pw <- function(mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, mu, kappa, gamma, pi){   #gamma as vector of offdiagonals
  mu.vonmises <- mu
  for(i in 1:4){
    mu <- cbind(mu1, mu2, mu3, mu4)
    assign(paste0("tmu", i), log(mu[,i]))    
  }
  for(i in 1:4){
    sigma <- cbind(sigma1, sigma2, sigma3, sigma4)
    assign(paste0("tsigma", i), log(sigma[,i]))    
  } 
  tmu <-  kappa*cos(mu.vonmises) # https://github.com/TheoMichelot/moveHMM/blob/master/R/n2w.R
  tkappa <- kappa*sin(mu.vonmises) 
  tgamma <- qlogis(gamma)
  tpi <- qlogis(pi)
  parvect <- c(tmu1, tmu2, tmu3, tmu4,tsigma1, tsigma2, tsigma3, tsigma4, tmu, tkappa, tgamma, tpi)
  return(parvect)
}



pw2pn <- function(parvect,N){     # parvect contains of mu, sigma as above, gamma
  mu1 <- exp(c(parvect[1:(1*N)]))
  mu2 <- exp(c(parvect[(N +1):(2*N)]))
  mu3 <- exp(c(parvect[(2*N +1):(3*N)]))
  mu4 <- exp(c(parvect[(3*N +1):(4*N)]))
  sigma1 <- exp(c(parvect[(1+4*N):(5*N)]))
  sigma2 <- exp(c(parvect[(1+5*N):(6*N)]))
  sigma3 <- exp(c(parvect[(1+6*N):(7*N)]))
  sigma4 <- exp(c(parvect[(1+7*N):(8*N)]))
  x <- parvect[(8*N+1):(9*N)]   # https://github.com/TheoMichelot/moveHMM/blob/master/R/w2n.R
  y <- parvect[(9*N+1):(10*N)]   
  mu <- Arg(x+1i*y)
  kappa <- sqrt(x^2+y^2)
  gamma <- diag(N) 
  if(N > 1){
    gamma[!gamma] <- exp(parvect[(10*N+1):(10*N+N*(N-1))])
    gamma <- gamma/apply(gamma,1,sum)           
  }
  delta <- solve(t(diag(N)-gamma+1),rep(1,N))
  pi <- plogis(parvect[((10*N+N*(N-1))+1):length(parvect)])
  return(list(mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, sigma1=sigma1, sigma2=sigma2, sigma3=sigma3, sigma4=sigma4,
              mu=mu, kappa=kappa, gamma=gamma, delta=delta, pi=pi))
}

#Next I will write my own function to fit an HMM based on roland code

L<-function(parvect, obs, N){ # parvect as working parameters with (divetim,maxdep,postdive,step,turn,matrix vectors for the off diagonals) variables with
  para <- pw2pn(parvect, N)
  mu.divetim <- para$mu1
  sigma.divetim <- para$sigma1
  mu.maxdep <- para$mu2
  sigma.maxdep <- para$sigma2
  mu.postdive.dur <- para$mu3
  sigma.postdive.dur <- para$sigma3
  mu.step <- para$mu4
  sigma.step <- para$sigma4
  mu.angle <- para$mu
  kappa <- para$kappa
  gamma <- matrix(para$gamma,N)
  delta <- c(para$delta)
  pi <- para$pi
  
  allprobs <- matrix(1,dim(obs)[1],N)
  #ind<-which(!is.na(obs$angle))    # angle has the most NA, but we can also exclude this observations before.
  ma <- matrix(1,dim(obs)[1],N) # including the zero inflation
  for (j in 1:N) {
      ind<-which(obs$postdive.dur==0)
      ma[ind,j] <- pgamma(obs$postdive.dur[ind]+0.5,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])
      ma[-ind,j]<- pgamma(obs$postdive.dur[-ind]+0.5,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])-
                   pgamma(obs$postdive.dur[-ind]-0.5,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])
  }
  for (j in 1:N){
    allprobs[,j] <-
      dgamma(obs$divetim, shape = mu.divetim[j]^2/sigma.divetim[j]^2, scale = sigma.divetim[j]^2/mu.divetim[j])*
      dgamma(obs$maxdep, shape = mu.maxdep[j]^2/sigma.maxdep[j]^2, scale = sigma.maxdep[j]^2/mu.maxdep[j])*
      ma[,j]*
      #dgamma(obs$postdive.dur, shape = mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2, scale = sigma.postdive.dur[j]^2/mu.postdive.dur[j])*
      dgamma(obs$step, shape = mu.step[j]^2/sigma.step[j]^2, scale = sigma.step[j]^2/mu.step[j])*
      dvm(obs$angle, mu.angle[j], kappa[j])
  }
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(obs$divetim)){
    foo <- phi%*%gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}


mle <- function(obs, mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, mu0, kappa0, gamma0, pi0, N){
  parvect <- pn2pw(mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, mu0, kappa0, gamma0, pi0)
  # obsl <- create_obslist(obs) 
  mod <- nlm(tsllk, parvect, obs, N, print.level = 2, iterlim = 1000, stepmax = 5) # replace L with tsllk
  pn <- pw2pn(mod$estimate, N)
  return(list(mu1 = pn$mu1, mu2 = pn$mu2, mu3 = pn$mu3, mu4 = pn$mu4, sigma1 = pn$sigma1, sigma2 = pn$sigma2, sigma3 = pn$sigma3, sigma4 = pn$sigma4,
              mu = pn$mu, kappa = pn$kappa, gamma = pn$gamma, delta = pn$delta, pi=pn$pi, mllk = mod$minimum))
}               


# time series likelihood
# a function which computes the likelihood for each time series after using the split function and adds the likelihoods
tsllk <- function(parvect,obs,N){
  # split the dataset to omit outliers in argument
  newdata <- splitset(k = 1800, argument = obs$postdive.dur, dataset = obs)
  # compute the likelihood for every subset and multiply them
  llk <- vector("list", length = length(newdata))
  llk[[1]] <- L(parvect, newdata[[1]], N)
  for (i in (2:length(newdata))){
    llk[[i]] <- llk[[(i-1)]] + L(parvect, newdata[[i]], N)
  }
  return(llk[[length(newdata)]])
}
  
  

  
#
fitmult <- function(obs, n_fits, N){
  modl <- list()
  for (i in 1:n_fits){
    modl[[i]] <- mle(obs, mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, mu0, kappa0, gamma0, N)
  }
  return(modl)
}              



################################################################################################################
### Testing and fitting
#mle(whaledata[(7:57), c(5, 7, 8, 13, 14)], c(20, 40), c(4, 8), c(15, 30), c(0.01, 0.03),
#    c(5, 10), c(1, 2), c(5, 20), c(0.05, 0.03), c(0.0), c(1,2), c(0.9, 0.8), 2)


#mle(whaledata[(600:750), c(5, 7, 8, 13, 14)], c(40, 100), c(10, 20), c(40, 65), c(0.01, 0.03),
#    c(10, 15), c(3, 5), c(15, 20), c(0.05, 0.03), c(0.0), c(1,2), c(0.9, 0.8), 2)
##

#obs <- nozeros[(7:2720), c(5, 7, 8, 13, 14)]
#obs <- whaledata[(7:2722), c(5, 7, 8, 13, 14)]
obs <- whaledata[(7:3583), c(5, 7, 8, 13, 14)]
obs <- obs[(obs$divetim > 1),]




#nozeros <- whaledivestats[(whaledivestats$postdive.dur > 0),]

###################
parvect <- pn2pw(c(20, 40), c(4, 8), c(15, 30), c(0.01, 0.03), c(5, 10), c(1, 2), c(5, 20), c(0.05, 0.03), c(0.0), c(1,2), c(0.9, 0.8))
para <- pw2pn(parvect, N)

parvect <- pn2pw(c(40, 100), c(10, 20), c(40, 65), c(0.01, 0.03), 
                 c(10, 15), c(3, 5), c(15, 20), c(0.05, 0.03), c(0.0), c(1,2), c(0.9, 0.8))
tsllk(parvect, obs, N)

N=2

###################
mle(obs, c(20, 40), c(4, 8), c(15, 30), c(0.01, 0.03),
    c(5, 10), c(1, 2), c(5, 20), c(0.05, 0.03), c(0.0), c(1,2), c(0.9, 0.8), c(0.2, 0.4), 2)
# andere Startwerte:
mod2 <- mle(obs, c(30, 200), c(20, 80), c(100, 600), c(10, 1000),
      c(25, 70), c(10, 20), c(50, 90), c(50, 500), c(0.0), c(1,2), c(0.9, 0.8), c(0.25,0.25), 2)

mle(obs, c(30, 200), c(20, 80), c(1, 6), c(0.01, 1),
    c(25, 70), c(10, 20), c(50, 90), c(0.05, 0.5), c(0.0), c(1,2), c(0.9, 0.8), c(0.25,0.25), 2)

mod3 <- mle(obs, c(30, 150, 250), c(20, 80, 110), c(50, 100, 600), c(10, 500, 1000),
    c(25, 60, 70), c(10, 20, 30), c(30, 90, 100), c(50, 100, 500), c(0.0), c(1, 2, 2.5), c(0.5, 0.4, 0.3, 0.2, 0.25, 0.1), c(0.25, 0.375, 0.375), 3)


#################
parvect<-c(3.205987e+00,  5.296187e+00 , 2.649179e+00,  4.358246e+00,
  4.214813e+00,  6.386805e+00, -4.230395e+00, -7.366967e-03,
  3.477124e+00 , 4.248472e+00,  2.386330e+00,  3.007819e+00,
  4.214810e+00,  4.504521e+00, -3.282264e+00, -6.886374e-01,
  1.133471e+00,  2.000357e+00, -3.202339e-04, -5.464697e-05,
  2.197357e+00,  1.179702e+00, -1.099844e+00, -1.098335e+00)

j=1
z<-seq(67.67101,67.69161,length=10000)
plot(z,dgamma(68, shape = mu.postdive.dur[j]^2/z^2, scale = z^2/mu.postdive.dur[j]))

z<-seq(0,1,length=10000)
#par(mfrow=(3,1))
plot(z,dgamma(z, shape =1, scale = 1))
plot(z,dgamma(z, shape =0.999, scale = 1))
plot(z,dgamma(z, shape =1.001, scale = 1),type="l")

dgamma(0, shape =1, scale = 1)


#####################################
mu01 <- c(30, 150, 250)
mu02 <- c(20, 80, 110)
mu03 <- c(50, 100, 600)
mu04 <- c(0.01, 0.5, 1)
sigma01 <- c(25, 60, 70)
sigma02 <- c(10, 20, 30)
sigma03 <- c(30, 90, 100)
sigma04 <- c(0.05, 0.1, 0.5)
mu.vm <- c(0.0)
kappa <- c(1, 2, 3)
gamma <- c(0.5, 0.4, 0.3)
pi0 <- c(0.25,0.25)
N <- 3
  
mle(obs, mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, mu.vm, kappa, gamma, pi0, N)



###################################################
## debugging ##
head(obs)
obs[,1]
hist(obs[,1])
hist(obs[,1],breaks=1000)
head(obs)
obs[which(obs[,1]==1),2]
obs[which(obs[,1]==1),3]
obs[which(obs[,1]==1),4]
hist(obs[which(obs[,1]==1),4])
max(obs[which(obs[,1]==1),4])
plot(obs[,4])
plot(obs[,4],type="h")

summary(obs[,1])
##################################################



## Plotting the 2-state-model solutions
hist(obs$divetim,probability=TRUE,breaks = 30)
z<-seq(0,600,by=0.01)
lines(z,mod2$delta[1]*dgamma(z,shape=mod2$mu1[1]^2/mod2$sigma1[1]^2,scale=mod2$sigma1[1]^2/mod2$mu1[1]),col='blue',lwd=2)
lines(z,mod2$delta[2]*dgamma(z,shape=mod2$mu1[2]^2/mod2$sigma1[2]^2,scale=mod2$sigma1[2]^2/mod2$mu1[2]),col='green',lwd=2)

hist(obs$maxdep,probability=TRUE,breaks = 30)
z<-seq(0,120,by=0.01)
lines(z,mod2$delta[1]*dgamma(z,shape=mod2$mu2[1]^2/mod2$sigma2[1]^2,scale=mod2$sigma2[1]^2/mod2$mu2[1]),col='blue',lwd=2)
lines(z,mod2$delta[2]*dgamma(z,shape=mod2$mu2[2]^2/mod2$sigma2[2]^2,scale=mod2$sigma2[2]^2/mod2$mu2[2]),col='green',lwd=2)

hist(obs$postdive.dur,probability=TRUE,breaks = 30)
z<-seq(0,5000,by=0.01)
lines(z,mod2$delta[1]*dgamma(z,shape=mod2$mu3[1]^2/mod2$sigma3[1]^2,scale=mod2$sigma3[1]^2/mod2$mu3[1]),col='blue',lwd=2)
lines(z,mod2$delta[2]*dgamma(z,shape=mod2$mu3[2]^2/mod2$sigma3[2]^2,scale=mod2$sigma3[2]^2/mod2$mu3[2]),col='green',lwd=2)

hist(obs$step,probability=TRUE,breaks = 30)
z<-seq(0,5000,by=0.01)
lines(z,mod2$delta[1]*dgamma(z,shape=mod2$mu4[1]^2/mod2$sigma4[1]^2,scale=mod2$sigma4[1]^2/mod2$mu4[1]),col='blue',lwd=2)
lines(z,mod2$delta[2]*dgamma(z,shape=mod2$mu4[2]^2/mod2$sigma4[2]^2,scale=mod2$sigma4[2]^2/mod2$mu4[2]),col='green',lwd=2)


## Plotting the 3-state-model solutions
hist(obs$divetim,probability=TRUE,breaks = 30)
z<-seq(0,600,by=0.01)
lines(z,mod3$delta[1]*dgamma(z,shape=mod3$mu1[1]^2/mod3$sigma1[1]^2,scale=mod3$sigma1[1]^2/mod3$mu1[1]),col='blue',lwd=2)
lines(z,mod3$delta[2]*dgamma(z,shape=mod3$mu1[2]^2/mod3$sigma1[2]^2,scale=mod3$sigma1[2]^2/mod3$mu1[2]),col='green',lwd=2)
lines(z,mod3$delta[3]*dgamma(z,shape=mod3$mu1[3]^2/mod3$sigma1[3]^2,scale=mod3$sigma1[3]^2/mod3$mu1[3]),col='red',lwd=2)

hist(obs$maxdep,probability=TRUE,breaks = 30)
z<-seq(0,120,by=0.01)
lines(z,mod3$delta[1]*dgamma(z,shape=mod3$mu2[1]^2/mod3$sigma2[1]^2,scale=mod3$sigma2[1]^2/mod3$mu2[1]),col='blue',lwd=2)
lines(z,mod3$delta[2]*dgamma(z,shape=mod3$mu2[2]^2/mod3$sigma2[2]^2,scale=mod3$sigma2[2]^2/mod3$mu2[2]),col='green',lwd=2)
lines(z,mod3$delta[3]*dgamma(z,shape=mod3$mu2[3]^2/mod3$sigma2[3]^2,scale=mod3$sigma2[3]^2/mod3$mu2[3]),col='red',lwd=2)

hist(obs$postdive.dur,probability=TRUE,breaks = 200, xlim=c(0,1000))
z<-seq(0,1000,by=0.01)
lines(z,mod3$delta[1]*dgamma(z,shape=mod3$mu3[1]^2/mod3$sigma3[1]^2,scale=mod3$sigma3[1]^2/mod3$mu3[1]),col='blue',lwd=2)
lines(z,mod3$delta[2]*dgamma(z,shape=mod3$mu3[2]^2/mod3$sigma3[2]^2,scale=mod3$sigma3[2]^2/mod3$mu3[2]),col='green',lwd=2)
lines(z,mod3$delta[3]*dgamma(z,shape=mod3$mu3[3]^2/mod3$sigma3[3]^2,scale=mod3$sigma3[3]^2/mod3$mu3[3]),col='red',lwd=2)

hist(obs$step,probability=TRUE,breaks = 30)
z<-seq(0,5000,by=0.01)
lines(z,mod3$delta[1]*dgamma(z,shape=mod3$mu4[1]^2/mod3$sigma4[1]^2,scale=mod3$sigma4[1]^2/mod3$mu4[1]),col='blue',lwd=2)
lines(z,mod3$delta[2]*dgamma(z,shape=mod3$mu4[2]^2/mod3$sigma4[2]^2,scale=mod3$sigma4[2]^2/mod3$mu4[2]),col='green',lwd=2)
lines(z,mod3$delta[3]*dgamma(z,shape=mod3$mu4[3]^2/mod3$sigma4[3]^2,scale=mod3$sigma4[3]^2/mod3$mu4[3]),col='red',lwd=2)

