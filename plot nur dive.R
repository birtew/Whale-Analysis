########################
## Creating the plots ##
########################

obs2 <- obs[(obs$postdive.dur < 1800),]
obs_new2 <- obs_new[(obs_new$postdive.dur < 1800),]


colpal <- c("#80cdc1", "gold", "#003c30") 




## Plotting the 2-state-model solutions
hist(obs_new2$divetim,probability=TRUE,breaks = 30, xlab = 'Divetime', main = 'State-dependent Distribution Divetime')
z<-seq(0,600,by=0.01)
lines(z,model2$delta[1]*dgamma(z,shape=model2$mu1[1]^2/model2$sigma1[1]^2,scale=model2$sigma1[1]^2/model2$mu1[1]),col=colpal[1],lwd=2)
lines(z,model2$delta[2]*dgamma(z,shape=model2$mu1[2]^2/model2$sigma1[2]^2,scale=model2$sigma1[2]^2/model2$mu1[2]),col=colpal[2],lwd=2)

hist(obs_new2$maxdep+5,probability=TRUE,breaks = 50, xlim = c(0, 120), xlab = 'Maximum depth', main = 'State-dependent Distribution Maximum Depth')
z<-seq(5,120,by=0.01)
lines(z,model2$delta[1]*dgamma(z-5,shape=model2$mu2[1]^2/model2$sigma2[1]^2,scale=model2$sigma2[1]^2/model2$mu2[1]),col=colpal[1],lwd=2)
lines(z,model2$delta[2]*dgamma(z-5,shape=model2$mu2[2]^2/model2$sigma2[2]^2,scale=model2$sigma2[2]^2/model2$mu2[2]),col=colpal[2],lwd=2)

hist(obs_new2$postdive.dur,probability=TRUE,breaks = 1000, xlim = c(0, 250), ylim = c(0, 0.05), xlab = 'Postdive Duration', main = 'State-dependent Distribution Postdive Duration')
z<-seq(0,250,by=0.01)
lines(z,model2$delta[1]*dgamma(z,shape=model2$mu3[1]^2/model2$sigma3[1]^2,scale=model2$sigma3[1]^2/model2$mu3[1]),col=colpal[1],lwd=2)
lines(z,model2$delta[2]*dgamma(z,shape=model2$mu3[2]^2/model2$sigma3[2]^2,scale=model2$sigma3[2]^2/model2$mu3[2]),col=colpal[2],lwd=2)



## Plotting the 3-state-model solutions
hist(obs_new2$divetim,probability=TRUE,breaks = 100, xlab = 'Divetime', main = 'State-dependent Distribution Divetime')
z<-seq(0,600,by=0.01)
lines(z,model3$delta[1]*dgamma(z,shape=model3$mu1[1]^2/model3$sigma1[1]^2,scale=model3$sigma1[1]^2/model3$mu1[1]),col=colpal[1],lwd=2)
lines(z,model3$delta[2]*dgamma(z,shape=model3$mu1[2]^2/model3$sigma1[2]^2,scale=model3$sigma1[2]^2/model3$mu1[2]),col=colpal[2],lwd=2)
lines(z,model3$delta[3]*dgamma(z,shape=model3$mu1[3]^2/model3$sigma1[3]^2,scale=model3$sigma1[3]^2/model3$mu1[3]),col=colpal[3],lwd=2)

# maxdep um 5 verschoben
hist(obs_new2$maxdep+5,probability=TRUE,breaks = 100, xlim = c(0, 120), xlab = 'Maximum depth', main = 'State-dependent Distribution Maximum Depth')
z<-seq(5,120,by=0.01)
lines(z,model3$delta[1]*dgamma(z-5,shape=model3$mu2[1]^2/model3$sigma2[1]^2,scale=model3$sigma2[1]^2/model3$mu2[1]),col=colpal[1],lwd=2)
lines(z,model3$delta[2]*dgamma(z-5,shape=model3$mu2[2]^2/model3$sigma2[2]^2,scale=model3$sigma2[2]^2/model3$mu2[2]),col=colpal[2],lwd=2)
lines(z,model3$delta[3]*dgamma(z-5,shape=model3$mu2[3]^2/model3$sigma2[3]^2,scale=model3$sigma2[3]^2/model3$mu2[3]),col=colpal[3],lwd=2)

hist(obs_new2$postdive.dur,probability=TRUE,breaks = 1000, xlim=c(0,250), ylim = c(0, 0.05), xlab = 'Postdive Duration', main = 'State-dependent Distribution Postdive Duration')
z<-seq(0,250,by=0.01)
lines(z,model3$delta[1]*dgamma(z,shape=model3$mu3[1]^2/model3$sigma3[1]^2,scale=model3$sigma3[1]^2/model3$mu3[1]),col=colpal[1],lwd=2)
lines(z,model3$delta[2]*dgamma(z,shape=model3$mu3[2]^2/model3$sigma3[2]^2,scale=model3$sigma3[2]^2/model3$mu3[2]),col=colpal[2],lwd=2)
lines(z,model3$delta[3]*dgamma(z,shape=model3$mu3[3]^2/model3$sigma3[3]^2,scale=model3$sigma3[3]^2/model3$mu3[3]),col=colpal[3],lwd=2)





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
  
  xi <- matrix(0,n,N)
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

states2 <- viterbi(obs_new2, c(model2$mu1), c(model2$mu2), c(model2$mu3), 
        c(model2$sigma1), c(model2$sigma2), c(model2$sigma3),
        c(model2$gamma), c(model2$delta), 2)



# Viterbi algorithm
# N = 3
states3 <- viterbi(obs_new2, c(model3$mu1), c(model3$mu2), c(model3$mu3), 
                   c(model3$sigma1), c(model3$sigma2), c(model3$sigma3),
                   c(model3$gamma), c(model3$delta), 3)


###
# Viterbi decoded plots
# neuer Vektor mit Farben je nach State
col2 <- states2
col2[which(col2 == 1)] <- colpal[1]
col2[which(col2 == 2)] <- colpal[2]

## 2 states
plot(obs_new2$divetim, type = 'h', col= col2, xlab = 'Dives', ylab = 'Divetime', main = 'Viterbi color-coded Divetimes')
plot(obs_new2$maxdep, type = 'h', col=col2)
plot(obs_new2$postdive.dur, type = 'h', col=col2)

# 3 states
col3 <- states3
col3[which(col3 == 1)] <- colpal[1]
col3[which(col3 == 2)] <- colpal[2]
col3[which(col3 == 3)] <- colpal[3]

plot(obs_new2$divetim, type = 'h', col= col3, xlab = 'Dives', ylab = 'Divetime', main = 'Viterbi color-coded Divetimes')
plot(obs_new2$maxdep, type = 'h', col=col3)
plot(obs_new2$postdive.dur, type = 'h', col=col3)



###







