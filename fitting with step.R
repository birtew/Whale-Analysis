require(CircStats) # for von Mises distribution
require(boot) # for logit
require(MASS)
require(reshape2)
require(stats)
require(moveHMM)


setwd("C:/Users/Patrick/Desktop/Whale")
#setwd("~/Uni/(M.Sc.) 3. Semester/Statistical Consulting/Minke whale project")
#whaledata <- read.csv("whale_data_cleaned.csv")
#whaledata <- read.csv("whale_data_cleaned_k2000.csv")
#whaledata <- read.csv("whale_data_cleaned_k2000_speed.csv")
whaledata <- read.csv("whale_data_cleaned_step.csv")

#obs <- whaledata[(7:3583), c(5, 7, 8, 13, 14)]
#obs <- obs[(obs$divetim > 1),]
obs <- whaledata[(7:3583), c(5, 7, 8, 15)]


## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
# mu & kappa: von Mises distr.
pn2pw <- function(mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, gamma, pi){   #gamma as vector of offdiagonals
  for(i in 1:4){
    mu <- cbind(mu1, mu2, mu3, mu4)
    assign(paste0("tmu", i), log(mu[,i]))    
  }
  for(i in 1:4){
    sigma <- cbind(sigma1, sigma2, sigma3, sigma4)
    assign(paste0("tsigma", i), log(sigma[,i]))    
  } 
  tgamma <- qlogis(gamma)
  tpi <- qlogis(pi)
  parvect <- c(tmu1, tmu2, tmu3, tmu4,tsigma1, tsigma2, tsigma3, tsigma4, tgamma, tpi)
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
  gamma <- diag(N) 
  if(N > 1){
    gamma[!gamma] <- exp(parvect[(8*N+1):(8*N+N*(N-1))])
    gamma <- gamma/apply(gamma,1,sum)           
  }
  delta <- solve(t(diag(N)-gamma+1),rep(1,N))
  pi <- plogis(parvect[((8*N+N*(N-1))+1):length(parvect)])
  return(list(mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, sigma1=sigma1, sigma2=sigma2, sigma3=sigma3, sigma4=sigma4,
              gamma=gamma, delta=delta, pi=pi))
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
  gamma <- matrix(para$gamma,N)
  delta <- c(para$delta)
  pi <- para$pi
  
  allprobs <- matrix(1,dim(obs)[1],N)
  #ind<-which(!is.na(obs$angle))    # angle has the most NA, but we can also exclude this observations before.
  ma <- matrix(1,dim(obs)[1],N) 
  for (j in 1:N) {
    ind<-which(obs$postdive.dur==0)
    ma[ind,j] <- pgamma(obs$postdive.dur[ind]+0.5,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])
    ma[-ind,j]<- pgamma(obs$postdive.dur[-ind]+0.5,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])-
      pgamma(obs$postdive.dur[-ind]-0.5,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])
  }
  ma2 <- matrix(1,dim(obs)[1],N) # including the zero inflation
  for (j in 1:N) {
    ind<-which(obs$step==0)
    ma2[ind,j] <- pi[j]
    ma2[-ind,j]<- (1-pi[j])*dgamma(obs$step[-ind], shape = mu.step[j]^2/sigma.step[j]^2, scale = sigma.step[j]^2/mu.step[j])
  }
  for (j in 1:N){
    allprobs[,j] <-
      (pgamma(obs$divetim +0.5, shape = mu.divetim[j]^2/sigma.divetim[j]^2, scale = sigma.divetim[j]^2/mu.divetim[j])-
         pgamma(obs$divetim -0.5, shape = mu.divetim[j]^2/sigma.divetim[j]^2, scale = sigma.divetim[j]^2/mu.divetim[j]))*
      dgamma(obs$maxdep, shape = mu.maxdep[j]^2/sigma.maxdep[j]^2, scale = sigma.maxdep[j]^2/mu.maxdep[j])*
      ma[,j]*
      ma2[,j]
    #dgamma(obs$postdive.dur, shape = mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2, scale = sigma.postdive.dur[j]^2/mu.postdive.dur[j])*
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



mle <- function(obs, mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, gamma0, pi0, N){
  parvect <- pn2pw(mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, gamma0, pi0)
  # obsl <- create_obslist(obs) 
  mod <- nlm(tsllk, parvect, obs, N, print.level = 2, iterlim = 1000, stepmax = 5) # replace L with tsllk
  pn <- pw2pn(mod$estimate, N)
  return(list(mu1 = pn$mu1, mu2 = pn$mu2, mu3 = pn$mu3, mu4 = pn$mu4, sigma1 = pn$sigma1, sigma2 = pn$sigma2, sigma3 = pn$sigma3, sigma4 = pn$sigma4,
              gamma = pn$gamma, delta = pn$delta, pi=pn$pi, mllk = mod$minimum))
}               


# ohne splitfunktion
mle <- function(obs, mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, gamma0, pi0, N){
  parvect <- pn2pw(mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, gamma0, pi0)
  # obsl <- create_obslist(obs) 
  mod <- nlm(L, parvect, obs, N, print.level = 2, iterlim = 1000, stepmax = 5) # replace L with tsllk
  pn <- pw2pn(mod$estimate, N)
  return(list(mu1 = pn$mu1, mu2 = pn$mu2, mu3 = pn$mu3, mu4 = pn$mu4, sigma1 = pn$sigma1, sigma2 = pn$sigma2, sigma3 = pn$sigma3, sigma4 = pn$sigma4,
              gamma = pn$gamma, delta = pn$delta, pi=pn$pi, mllk = mod$minimum))
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



# Splitset function
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



#
fitmult <- function(obs, n_fits, N){
  modl <- list()
  for (i in 1:n_fits){
    modl[[i]] <- mle(obs, mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, mu0, kappa0, gamma0, N)
  }
  return(modl)
}              



################################################################################################################
### Testing and fitting ###
###########################
# 2 states
whale_mod2withstep <- mle(obs, c(30, 200), c(20, 80), c(100, 600), c(10, 1000),
                          c(25, 70), c(10, 20), c(50, 90), c(50, 500), c(0.9, 0.8), c(0.25, 0.25), 2)

# 3 states
whale_mod3withstep <- mle(obs, c(30, 150, 250), c(20, 80, 110), c(50, 100, 600), c(10, 500, 1000),
                          c(25, 60, 70), c(10, 20, 30), c(30, 90, 100), c(50, 100, 500), c(0.5, 0.4, 0.3, 0.2, 0.25, 0.1), c(0.25, 0.25, 0.25), 3)


## zufällig gezogene Startwerte:
# 2 states
whale_rmod2withstep <- mle(obs, c(runif(2, 1, 250)), c(runif(2, 1, 120)), c(runif(2, 1, 600)), c(runif(2, 0.01, 100)),
                           c(runif(2, 1, 100)), c(runif(2, 1, 50)), c(runif(2, 1, 100)), c(runif(2, 1, 500)), 
                           c(runif(2, 0, 1)), c(runif(2, 0, 1)), 2)

# 3 states
whale_rmod3withstep <- mle(obs, c(runif(3, 1, 250)), c(runif(3, 1, 120)), c(runif(3, 1, 600)), c(runif(3, 0.01, 100)),
                           c(runif(3, 1, 100)), c(runif(3, 1, 50)), c(runif(3, 1, 100)), c(runif(3, 1, 500)), 
                           c(runif(6, 0, 1)), c(runif(3, 0, 1)), 3)


## Dealing with lokal maxima
# 2 states:
llks<-rep(NA,10)
mods2_step <- vector("list")
for (k in 1:10){
  mods2_step[[k]] <- mle(obs, c(runif(2, 1, 250)), c(runif(2, 1, 120)), c(runif(2, 1, 600)), c(runif(2, 0.01, 100)),
                   c(runif(2, 1, 100)), c(runif(2, 1, 50)), c(runif(2, 1, 100)), c(runif(2, 1, 500)), 
                   c(runif(2, 0, 1)), c(runif(2, 0, 1)), 2)
  llks[k] <- -mods2_step[[k]]$mllk #minimum
}
model2_step <- mods2_step[[which.max(llks)]]

# 3 states:
llks<-rep(NA,10)
mods3_step <- vector("list")
for (k in 1:10){
  mods3_step[[k]] <- mle(obs, c(runif(3, 1, 250)), c(runif(3, 1, 120)), c(runif(3, 1, 600)), c(runif(3, 0.01, 100)),
                         c(runif(3, 1, 100)), c(runif(3, 1, 50)), c(runif(3, 1, 100)), c(runif(3, 1, 500)), 
                         c(runif(6, 0, 1)), c(runif(3, 0, 1)), 3)
  llks[k] <- -mods3_step[[k]]$mllk #minimum
}
model3_step <- mods3_step[[which.max(llks)]]


#####################################
mu01 <- c(30, 150, 250)
mu02 <- c(20, 80, 110)
mu03 <- c(50, 100, 600)
mu04 <- c(0.01, 0.5, 1)
sigma01 <- c(25, 60, 70)
sigma02 <- c(10, 20, 30)
sigma03 <- c(30, 90, 100)
sigma04 <- c(50, 100, 500)
mu.vm <- c(0.0)
kappa <- c(1, 2, 3)
gamma <- c(0.5, 0.4, 0.3)
pi0 <- c(0.25,0.25)
N <- 3

mle(obs, mu01, mu02, mu03, mu04, sigma01, sigma02, sigma03, sigma04, mu.vm, kappa, gamma, pi0, N)





###########################################################################
# Debugging
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


#########
speed <- obs$step / obs$divetim
boxplot(speed)
which(speed > 11)
summary(speed)
speed_kmh <- speed * 3.6
summary(speed_kmh)

