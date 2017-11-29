require(CircStats) # for von Mises distribution
require(boot) # for logit
require(MASS)
require(reshape2)
require(stats)
require(moveHMM)


setwd("C:/Users/Patrick/Desktop/Whale")
whaledata <- read.csv("whale_data_cleaned.csv")



 
## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
# mu & kappa: von Mises distr.
pn2pw <- function(mu1, mu2, mu3, mu4, sigma1, sigma2, sigma3, sigma4, mu, kappa, gamma){   #gamma as vector of offdiagonls
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
  parvect <- c(tmu1, tmu2, tmu3, tmu4,tsigma1, tsigma2, tsigma3, tsigma4,tmu,tkappa,tgamma)
  return(parvect)
}



pw2pn <- function(parvect,N){     # parvct contains of mu, sigma as above, gamma
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
    gamma[!gamma] <- exp(parvect[(10*N+1):(length(parvect))])
    gamma <- gamma/apply(gamma,1,sum)           
  }
  delta <- solve(t(diag(N)-gamma+1),rep(1,N))
  return(list(mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,sigma4=sigma4,
              mu=mu,kappa=kappa,gamma=gamma,delta=delta))
}

#Next I will write my own function to fit an HMM based on roland code

L<-function(parvect,obs,N){ # parvect as working parameters with (divetim,maxdep,postdive,step,turn,matrix vectors for the off diagonals) variables with
  para <- pw2pn(parvect, N)
  mu.divetim <-para$mu1
  sigma.divetim <- para$sigma1
  mu.maxdep <-para$mu2
  sigma.maxdep <- para$sigma2
  mu.postdive.dur <- para$mu3
  sigma.postdive.dur <- para$sigma3
  mu.step <- para$mu4
  sigma.step <- para$sigma4
  mu.angle <- para$mu
  kappa <- para$kappa
  gamma <- matrix(para$gamma,N)
  delta <- c(para$delta)
  
  
  allprobs <- matrix(1,dim(obs)[1],N)
  #ind<-which(!is.na(obs$angle))    #angle has the most NA, but we can also exclude this observations before.
  for (j in 1:N){
    allprobs[,j] <-
      dgamma(obs$divetim,shape=mu.divetim[j]^2/sigma.divetim[j]^2,scale=sigma.divetim[j]^2/mu.divetim[j])*
      dgamma(obs$maxdep,shape=mu.maxdep[j]^2/sigma.maxdep[j]^2,scale=sigma.maxdep[j]^2/mu.maxdep[j])*
      dgamma(obs$postdive.dur,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])*
      dgamma(obs$step,shape=mu.step[j]^2/sigma.step[j]^2,scale=sigma.step[j]^2/mu.step[j])*
      dvm(obs$angle,mu.angle[j],kappa[j])
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


mle <- function(obs,mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,mu0,kappa0,gamma0,N){
  parvect <- pn2pw(mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,mu0,kappa0,gamma0)
  # obsl <- create_obslist(obs) 
  mod <- nlm(L,parvect,obs,N,print.level=2,iterlim=1000,stepmax=5)
  pn <- pw2pn(mod$estimate,N)
  return(list(mu1=pn$mu1,mu2=pn$mu2,mu3=pn$mu3,mu4=pn$mu4,sigma1=pn$sigma1,sigma2=pn$sigma2,sigma3=pn$sigma3,sigma4=pn$sigma4,
              mu=pn$mu,kappa=pn$kappa,gamma=pn$gamma,delta=pn$delta,mllk=mod$minimum))
}               



fitmult <- function(obs,n_fits,N){
  modl <- list()
  for (i in 1:n_fits){
    modl[[i]] <- mle(obs,mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,mu0,kappa0,gamma0,N)
  }
  return(modl)
}              


###
mle(whaledata[(7:57),c(5,7,8,13,14)], c(20,40), c(4, 8), c(15,30), c(0.01, 0.03),
    c(5, 10), c(1,2), c(5, 20), c(0.05, 0.03), c(0.0), c(1,2), c(0.9, 0.8), 2)


