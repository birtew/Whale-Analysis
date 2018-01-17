############################################
# Pseudoresiduen für divetime und maxdepth #

lforward <- function(x,mu,sigma,Gamma,delta,N){
  n <- length(x)
  lalpha <- matrix (NA,N,n)
  shape <- mu^2/sigma^2
  scale <- sigma^2/mu
  allprobs<-matrix(1,n,N)
  ind<-which(!is.na(x))
  for (j in 1:N){
    allprobs[ind,j]<-dgamma(x[ind],shape=shape[j],scale=scale[j])
  }
  foo <- delta*allprobs[1,]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[,1] <- lscale+log(foo)
  for (i in 2:n)
  {
    foo <- foo%*%Gamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo)+lscale
  }
  return(lalpha)
}


PseudoRes<- function(x,mu,sigma,Gamma,delta,N){
  n <-length(x)
  la <- lforward(x,mu,sigma,Gamma,delta,N)
  shape <- mu^2/sigma^2
  scale <- sigma^2/mu
  Res<-rep(NA,n)
  P<-matrix(NA,n,N)
  for (j in 1:N){
    P[,j]<-pgamma(x,shape=shape[j],scale=scale[j])
  }
  Res[1]<-qnorm(t(delta)%*%P[1,])
  for (i in 2:n){
    c<-max(la[,i-1])
    a<-exp(la[,i-1]-c)
    Res[i]<-qnorm(t(a)%*%(Gamma/sum(a))%*%P[i,])
  }
  return(Res)
}

# Pseudoresiduen für postdivedur. Achtung: Große Werte in Postdivedur müssen eliminiert werden, Beobachtungen müssen ohne NA sein, aber das ist bei uns gegeben

lforwardpostdive <- function(x,mu,sigma,Gamma,delta,N){
  n <- length(x)
  lalpha <- matrix (NA,N,n)
  shape <- mu^2/sigma^2
  scale <- sigma^2/mu
  allprobs<-matrix(1,n,N)
  ind<-which(x==0)
  for (j in 1:N) {
    allprobs[ind,j] <- pgamma(x[ind]+0.5,shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])
    allprobs[-ind,j]<- pgamma(x[-ind]+0.5,shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])-
      pgamma(x[-ind]-0.5,shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])
  }
  foo <- delta*allprobs[1,]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[,1] <- lscale+log(foo)
  for (i in 2:n)
  {
    foo <- foo%*%Gamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo)+lscale
  }
  return(lalpha)
}



PseudoRespostdive <- function(x,mu,sigma,Gamma,delta,N){
  n <-length(x)
  la <- lforwardpostdive(x,mu,sigma,Gamma,delta,N)
  shape <- mu^2/sigma^2
  scale <- sigma^2/mu
  Res<-rep(NA,n)
  P<-matrix(NA,n,N)
  ind<-which(x==0)
  for (j in 1:N){
    P[ind,j] <- pgamma(x[ind]+0.5,shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])
    P[-ind,j]<- pgamma(x[-ind]+0.5,shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])-
      pgamma(x[-ind]-0.5,shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])
  }
  Res[1]<-qnorm(t(delta)%*%P[1,])
  for (i in 2:n){
    c<-max(la[,i-1])
    a<-exp(la[,i-1]-c)
    Res[i]<-qnorm(t(a)%*%(Gamma/sum(a))%*%P[i,])
  }
  return(Res)
}

################################################################################################
##############################
## 2 States
# Residuen für divetim
Res_divetim <- PseudoRes(obs_new2$divetim, model2$mu1, model2$sigma1, matrix(c(model2$gamma), 2), model2$delta, 2)
hist(Res_divetim)
acf(Res_divetim,na.action=na.pass)
qqnorm(Res_divetim)
qqline(Res_divetim, col = 2)

# Residuen für maxdep
Res_maxdep <- PseudoRes(obs_new2$maxdep, model2$mu2, model2$sigma2, matrix(c(model2$gamma), 2), model2$delta, 2)
hist(Res_maxdep)
acf(Res_maxdep,na.action=na.pass)
qqnorm(Res_maxdep)
qqline(Res_maxdep, col = 2)

# Residuen für postdive.dur
Res_postdive <- PseudoRespostdive(obs_new2$postdive.dur, model2$mu3, model2$sigma3, matrix(c(model2$gamma), 2), model2$delta, 2)
hist(Res_postdive)
acf(Res_postdive,na.action=na.pass)
qqnorm(Res_postdive)
qqline(Res_postdive, col = 2)


##############################
## 3 States
# Residuen für divetim
Res_divetim <- PseudoRes(obs_new2$divetim, model3$mu1, model3$sigma1, matrix(c(model3$gamma), 3), model3$delta, 3)
hist(Res_divetim)
acf(Res_divetim,na.action=na.pass)
qqnorm(Res_divetim)
qqline(Res_divetim, col = 2)

# Residuen für maxdep
Res_maxdep <- PseudoRes(obs_new2$maxdep, model3$mu2, model3$sigma2, matrix(c(model3$gamma), 3), model3$delta, 3)
hist(Res_maxdep)
acf(Res_maxdep,na.action=na.pass)
qqnorm(Res_maxdep)
qqline(Res_maxdep, col = 2)

# Residuen für postdive.dur
Res_postdive <- PseudoRespostdive(obs_new2$postdive.dur, model3$mu3, model3$sigma3, matrix(c(model3$gamma), 3), model3$delta, 3)
hist(Res_postdive)
acf(Res_postdive,na.action=na.pass)
qqnorm(Res_postdive)
qqline(Res_postdive, col = 2)


