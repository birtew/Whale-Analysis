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
pn2pw <- function(mu1,sigma1,mu2,sigma2,mu3,sigma3,mu4,sigma4,mu,kappa,gamma,N){   
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
  tgamma <- NULL
  if(N>1){
    foo <- log(gamma/diag(gamma))
    tgamma <- as.vector(foo[!diag(N)])
  }
  parvect <- c(tmu1,tsigma1,tmu2,tsigma2,tmu3,tsigma3,tmu4,tsigma4,tmu,tkappa,tgamma)
  return(parvect)
}

pw2pn <- function(parvect,N){
  mu <- exp(c(parvect[1:N],parvect[2*N+1:3*N],parvect[4*N+1:3*N],parvect[6*N+1:7*N])
  for(i in 1:4) {
    assign(paste0("mu", i), mu[(i*N-N+1):(i*N)])
  }   
  sigma <- exp(c(parvect[1+N:2*N],parvect[3*N+1:4*N],parvect[5*N+1:6*N],parvect[7*N+1:8*N])
  for(i in 1:4) {
    assign(paste0("sigma", i), sigma[(i*N-N+1):(i*N)])
  } 
  x <- parvect[(8*N+1):(9*N)]   # https://github.com/TheoMichelot/moveHMM/blob/master/R/w2n.R
  y <- parvect[(9*N+1):(10*N)]   
  mu <- Arg(x+1i*y)
  kappa <- sqrt(x^2+y^2)
  gamma <- diag(N)
  if(N>1){
    gamma[!gamma] <- exp(parvect[(10*N+1):(length(parvect))])
    gamma <- gamma/apply(gamma,1,sum)           
  }
  delta <- solve(t(diag(N)-gamma+1),rep(1,N))
  return(list(mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,sigma4=sigma4,
              mu=mu,kappa=kappa,gamma=gamma,delta=delta))
}

#Next I will write my own function to fit an HMM based on roland code

L<-function(theta.star,x,N){ # theta.star as working parameters with (divetim,maxdep,postdive,step,turn,matrix vectors for the off diagonals) variables with
  mu.divetim <- exp(theta.star[1:N])
  sigma.divetim <- exp(theta.star[N+1:2*N])
  mu.maxdep <- exp(theta.star[2*N+1:3*N])
  sigma.maxdep <- exp(theta.star[3*N+1:4*N])
  mu.postdive.dur <- exp(theta.star[4*N+1:5*N])
  sigma.postdive.dur <- exp(theta.star[5*N+1:6*N])
  mu.step <- exp(theta.star[6*N+1:7*N])
  sigma.step <- exp(theta.star[7*N+1:8*N])
  mu.turn <- kappa*cos(theta.star[8*N+1:9*N])   # https://github.com/TheoMichelot/moveHMM/blob/master/R/n2w.R
  sigma.turn <- kappa*sin(theta.star[9*N+1:10*N])
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[10*N+1:10*N+((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  allprobs <- matrix(1,dim(x)[1],N)
  ind<-which(!is.na(x$angle))    #angle has the most NA, but we can also exclude this observations before.
  for (j in 1:N){
    allprobs[,j] <-
      dgamma(x$divetim,shape=mu.divetim[j]^2/sigma.divetim[j]^2,scale=sigma.divetim[j]^2/mu.divetim[j])*
      dgamma(x$maxdep,shape=mu.maxdep[j]^2/sigma.maxdep[j]^2,scale=sigma.maxdep[j]^2/mu.maxdep[j])*
      dgamma(x$postdive.dur,shape=mu.postdive.dur[j]^2/sigma.postdive.dur[j]^2,scale=sigma.postdive.dur[j]^2/mu.postdive.dur[j])*
      dgamma(x$step,shape=mu.step[j]^2/sigma.step[j]^2,scale=sigma.step[j]^2/mu.step[j])*
      dvm(x$turn,mu.turn[j],kappa[j])
  }
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x$divetim)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}
