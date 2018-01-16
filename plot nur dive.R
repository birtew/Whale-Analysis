#######################
## Creating the plots

## Plotting the 2-state-model solutions
hist(obs$divetim,probability=TRUE,breaks = 30)
z<-seq(0,600,by=0.01)
lines(z,mod2$delta[1]*dgamma(z,shape=mod2$mu1[1]^2/mod2$sigma1[1]^2,scale=mod2$sigma1[1]^2/mod2$mu1[1]),col='blue',lwd=2)
lines(z,mod2$delta[2]*dgamma(z,shape=mod2$mu1[2]^2/mod2$sigma1[2]^2,scale=mod2$sigma1[2]^2/mod2$mu1[2]),col='green',lwd=2)

hist(obs$maxdep,probability=TRUE,breaks = 30)
z<-seq(0,120,by=0.01)
lines(z,mod2$delta[1]*dgamma(z,shape=mod2$mu2[1]^2/mod2$sigma2[1]^2,scale=mod2$sigma2[1]^2/mod2$mu2[1]),col='blue',lwd=2)
lines(z,mod2$delta[2]*dgamma(z,shape=mod2$mu2[2]^2/mod2$sigma2[2]^2,scale=mod2$sigma2[2]^2/mod2$mu2[2]),col='green',lwd=2)

hist(obs$postdive.dur,probability=TRUE,breaks = 100, xlim = c(0, 1500))
z<-seq(0,1500,by=0.01)
lines(z,mod2$delta[1]*dgamma(z,shape=mod2$mu3[1]^2/mod2$sigma3[1]^2,scale=mod2$sigma3[1]^2/mod2$mu3[1]),col='blue',lwd=2)
lines(z,mod2$delta[2]*dgamma(z,shape=mod2$mu3[2]^2/mod2$sigma3[2]^2,scale=mod2$sigma3[2]^2/mod2$mu3[2]),col='green',lwd=2)

hist(obs$step,probability=TRUE,breaks = 100, xlim = c(0, 1500))
z<-seq(0,1500,by=0.01)
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

hist(obs$step,probability=TRUE,breaks = 200, xlim = c(0, 1000))
z<-seq(0,1000,by=0.01)
lines(z,mod3$delta[1]*dgamma(z,shape=mod3$mu4[1]^2/mod3$sigma4[1]^2,scale=mod3$sigma4[1]^2/mod3$mu4[1]),col='blue',lwd=2)
lines(z,mod3$delta[2]*dgamma(z,shape=mod3$mu4[2]^2/mod3$sigma4[2]^2,scale=mod3$sigma4[2]^2/mod3$mu4[2]),col='green',lwd=2)
lines(z,mod3$delta[3]*dgamma(z,shape=mod3$mu4[3]^2/mod3$sigma4[3]^2,scale=mod3$sigma4[3]^2/mod3$mu4[3]),col='red',lwd=2)



# Viterbi algorithm
# N = 2
viterbi<-function(obs, mu1, mu2, mu3, sigma1, sigma2, sigma3, gamma, delta, N){
  n <- length(obs$divetim)
  
  mu.divetim <- mu1
  sigma.divetim <- sigma1
  mu.maxdep <- mu2
  sigma.maxdep <- sigma2
  mu.postdive.dur <- mu3
  sigma.postdive.dur <- sigma3
  gamma <- matrix(gamma,N)
  delta <- c(delta)
  
  allprobs <- matrix(1,dim(obs)[1],N)
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
      ma[,j]
  }
  
  xi <- matrix(0,n,2)
  foo <- delta*allprobs[1,]
  xi[1,] <- foo/sum(foo)
  for (i in 2:n){
    foo <- apply(xi[i-1,]*gamma,2,max)*allprobs[i,]
    xi[i,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for (i in (n-1):1){
    iv[i] <- which.max(gamma[,iv[i+1]]*xi[i,])
  }
  iv
}

obs2 <- obs[(obs$postdive.dur < 1800),]

states <- viterbi(obs2, unlist(whale_dive_rmod2$mu1), unlist(whale_dive_rmod2$mu2), unlist(whale_dive_rmod2$mu3), 
        unlist(whale_dive_rmod2$sigma1), unlist(whale_dive_rmod2$sigma2), unlist(whale_dive_rmod2$sigma3),
        unlist(whale_dive_rmod2$gamma), unlist(whale_dive_rmod2$delta), 2)

#######
plot(obs2$divetim[1:100])

#plot the results:
b <- seq(-bm,bm,length=m+1)                 
h <- b[2]-b[1]                              
bstar <- (b[-1]+b[-(m+1)])*0.5              
vola <- bstar[states] # Viterbi-decoded volatility levels

par(mfrow=c(2,1))
plot(obs2$divetim,type="l",xlab="time",ylab="daily returns",main="time series of Deutsch Bank share returns",cex.main=1.2)
plot(vola,type="l",xlab="time",ylab="decoded (log-)volatilities")






