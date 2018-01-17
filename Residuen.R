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
  la <- lforward(x,mu,sigma,Gamma,delta,N)
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



#example for N=2 
testres <- mle(obs, c(30, 200), c(20, 80), c(100, 600),
                       c(25, 70), c(10, 20), c(50, 90), c(0.9, 0.8), 2)
                       
x <- obs[which(obs$postdive.dur < 1800),]$postdive.dur
mu <- testres$mu3
sigma <- testres$sigma3
Gamma <- matrix(c(testres$gamma),2)
delta <- testres$delta
N=2

PseudoRespostdive(x,mu,sigma,Gamma,delta,N)

hist(Res)

acf(Res,na.action=na.pass)
qqnorm(Res)
qqline(Res, col = 2)

# Pseudoresiduen für divetime und masdepth. Achtung: Große Werte in Postdivedur müssen eliminiert werden


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




#example for N=2 
testres <- mle(obs, c(30, 200), c(20, 80), c(100, 600),
                       c(25, 70), c(10, 20), c(50, 90), c(0.9, 0.8), 2)
                       
x <- obs[which(obs$postdive.dur < 1800),]$divetim
mu <- testres$mu1
sigma <- testres$sigma1
Gamma <- matrix(c(testres$gamma),2)
delta <- testres$delta
N=2

PseudoRes(x,mu,sigma,Gamma,delta,N)

hist(Res)

acf(Res,na.action=na.pass)

qqnorm(Res)
qqline(Res, col = 2)
