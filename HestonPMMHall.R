library(latex2exp)
library(matrixcalc)
library(emdbook)
library(MASS)
library(latex2exp)
library(Metrics)
library(MLmetrics)
library(future.apply)
library(future)
library(gmodels)






set.seed(3005)  #2 was good! as was 3 1116 with 1mill results 'not bad'
parameters <- c(1,5,0.15,0.1,0.7) #mu kappa theta xi rho
Hdata <- hestonSim(hParams=parameters,S0=100,V0=0.05,T=1,deltat=1/(10*52))
Hprice <- Hdata[[1]][1+(0:51)*10] #thin to get data points every 1/52 time units
Hvols <- Hdata[[2]][1+(0:51)*10]
data <- matrix(c(Hprice, Hvols), nrow=2, byrow=TRUE)
plot(ts(Hprice,deltat=1/52,start=0))
plot(ts(Hvols,deltat=1/52,start=0))
test <- heston_pmmh(100000,data,0.1)


plot(ts(test[1,],start=0,deltat=1),ylab=TeX("\\mu"),xlab="Iteration",main=TeX("$\\mu$ Trace plot"))

plot(ts(exp(test[2,]),start=0,deltat=1),ylab=TeX("\\kappa"),xlab="Iteration",main=TeX("$\\kappa$ Trace plot"))


plot(ts(exp(test[3,]),start=0,deltat=1),ylim=c(0,0.2),ylab=TeX("\\theta"),xlab="Iteration",main=TeX("$\\theta$ Trace plot"))


plot(ts(exp(test[4,]),start=0,deltat=1),ylab=TeX("\\xi"),xlab="Iteration",main=TeX("$\\xi$ Trace plot"))


png('rhochain.png')
plot(ts(test[5,],start=0,deltat=1),ylab=TeX("\\rho"),xlab="Iteration",main=TeX("$\\rho$ Trace plot"))


hist(test[1,],freq=F,main=TeX("$\\mu$ Histogram"),xlab=TeX("\\mu"), axes=FALSE,breaks=100,ylim=c(0,3),mgp=c(1.5,1,0),xlim=c(0,1.5))
axis(1,at = seq(0,1.5,0.5),labels = TRUE,pos = 0)
axis(2,pos = 0)
mean(test[1,])
abline(v=mean(test[1,]),col="blue",lwd=4)
abline(v=parameters[1],col="red",lwd=4)

hist(exp(test[2,]),freq=F,main=TeX("$\\kappa$ Histogram"),xlab=TeX("\\kappa"), axes=FALSE,breaks=100,ylim=c(0,3),mgp=c(1.5,1,0),xlim=c(4,8))
exp(mean(test[2,]))
axis(1,at = seq(4,8,0.5),labels = TRUE,pos = 0)
axis(2,pos = 4)
abline(v=exp(mean(test[2,])),col="blue",lwd=4)
abline(v=parameters[2],col="red",lwd=4)

hist(exp(test[3,]),freq=F,main=TeX("$\\theta$ Histogram"),xlab=TeX("\\theta"),axes=FALSE, breaks=1000,ylim=c(0,25),mgp=c(1.5,1,0),xlim=c(0.1,0.2))
axis(1,at = seq(0.1,0.2,0.05),labels = TRUE,pos = 0)
axis(2,pos = 0.1)
exp(mean(test[3,]))
abline(v=exp(mean(test[3,])),col="blue",lwd=4)
abline(v=parameters[3],col="red",lwd=4)

hist(exp(test[4,]),freq=F,main=TeX("$\\xi$ Histogram"),xlab=TeX("\\xi"),axes=FALSE, breaks=100,ylim=c(0,25),mgp=c(1.5,1,0),xlim=c(0.075,0.15))
axis(1,at = seq(0.075,0.15,0.05),labels = TRUE,pos = 0)
axis(2,pos = 0.075)
exp(mean(test[4,]))
abline(v=exp(mean(test[4,])),col="blue",lwd=4)
abline(v=parameters[4],col="red",lwd=4)

hist(test[5,],freq=F,main=TeX("$\\rho$ Histogram"),xlab=TeX("\\rho"), breaks=100,ylim=c(0,10),mgp=c(1.5,1,0),xlim=c(0.5,1),axes=FALSE)
axis(1,at = seq(0.5,1,0.5),labels = TRUE,pos = 0)
axis(2,pos = 0.5)
mean(test[5,])
abline(v=mean(test[5,]),col="blue",lwd=4)
abline(v=parameters[5],col="red",lwd=4)


##Change index of test and main accordingly
plot(acf(test[2,],lag.max=1500),main=TeX("$\\kappa$  Corellogram"),ylab=TeX("$\\rho_k"))
abline(h=0.006,col="red")
abline(h=-0.006,col="red")
dataf= data.frame(test[1,],test[2,],test[3,],test[4,],test[5,])




heston_pmmh <- function(n, X,l){
  
  theta <- array(0,dim=c(5,n)) #mu kappa theta xi rho
  theta[,1][1] <- 0.1
  theta[,1][2] <- log(0.1)
  theta[,1][3] <- log(0.1)
  theta[,1][4] <- log(0.1)
  theta[,1][5] <- 0.1
  
  count <- 0
  input_params <- vector()
  input_params[1] <- theta[,1][1]
  input_params[2:4] <- exp(theta[,1][2:4])
  input_params[5] <- theta[,1][5]
  proposal_prev <- loglikelihood_heston(X,params=input_params) +log_prior_heston(theta[,1])
  for (i in 2:n){
    if(i%%100==0){
      print(i)
      #print(theta[,i-1])
    }
    proposal_theta <- rnorm(5,theta[,i-1],l)
    input_proposal <- vector()
    input_proposal[1] <- proposal_theta[1]
    input_proposal[2:4] <- exp(proposal_theta[2:4])
    input_proposal[5] <- proposal_theta[5]
    target_proposal <- loglikelihood_heston(X,params=input_proposal) +log_prior_heston(proposal_theta)
    alpha <- target_proposal-proposal_prev
    if(is.na(alpha)){
      print(i)
      print(proposal_theta)
      print(proposal_prev)
      print(target_proposal)
    }
    if (log(runif(1))<alpha){
      theta[,i] <- proposal_theta
      proposal_prev <- target_proposal
      count <- count + 1
    }
    else {
      theta[,i] <- theta[,i-1]
    }
  }
  print(count/(n-1))
  #return(count/(n-1))
  return(theta)
}


ResidualbridgejointdensityMV=function(x,T,eta,afun=alphaLV_sabr,bfun=betaLV_sabr,deltat,params)
{
  n = T/deltat
  ll=0
  full <- x
  x <- x-eta
  x1=x[n+1,]
  
  for (i in 1:(n-1)){
    t=i*deltat
    ll=ll+lgpdf(x[i+1,],((T-t)*x[i,]+x1*deltat)/(T-t+deltat),((T-t)/(T-t+deltat))*deltat*bfun(full[i,],params))
  }
  return(ll)
}

heston_residual <- function(x0,xT,eta,deltat,T,afun,bfun,params){ #bfun is just the heston mdb #afun is bit different
  n = T/deltat
  d=length(x0) #d is state dimension
  x=matrix(0,ncol=d,nrow=n+1) # columns correspond to each component
  rb <- matrix(0,ncol=d,nrow=n+1)
  x[1,]=x0
  x[n+1,]=xT
  nT <- eta[n+1,]
  for(i in 1:(n-1))
  {
    delta_eta <- (eta[i+1,]-eta[i,])/deltat   #The eta defined in thesis and 'Improved bridge constructs'
    tm = (i-1)*deltat
    t = i*deltat
    #print(x[i,])
    x[i+1,] = (x[i,]+ delta_eta*deltat +(xT-x[i,]-nT+eta[i,])*(deltat)/(T-tm)+sqrt((T-t)/(T-tm))*t(chol(bfun(x[i,],params)))%*%rnorm(2,0,sqrt(deltat)))
    # for(j in 1:d){
    #   if(x[i+1,j]<0.01){
    #     x[i+1,j]=0.01 #fudge, to avoid negative values
    #   }
    # }
  }
  return(x)
}

heston_det_func <- function(params,n0,T,deltat){ #eta deterministic function. params = mu kappa theta
  mu <- params[1]
  kappa <- params[2]
  theta <- params[3]
  eta <- array(0,dim=c((T/deltat)+1,2))
  eta_1 <- seq(0,T,by=deltat)
  eta_1 <- n0[1]*exp(mu*eta_1)
  eta_2 <- seq(0,T,by=deltat)
  eta_2 <- theta-exp(-kappa*eta_2)*(theta-n0[2])
  eta[,1] <- eta_1
  eta[,2] <- eta_2
  #print(eta[,2])
  return(eta)
}

wr_heston <- function(x0, xT, T=10, deltat=1/52, N=1, params,afun=alpha_heston,bfun=beta_heston,wr=FALSE) #set wr=TRUE to perform weighted resampling. Else importance weights are returned
{
  paths <- array(0,dim=c((T/deltat)+1,2,N))
  logweights <- rep(0,N)
  paths[1] <- 2
  for (i in 1:N)
  {
    heston_det <- heston_det_func(params,x0,T,deltat)
    heston_diffusion <- heston_residual(x0, xT,heston_det, deltat, T,afun,bfun,params)
    path <- heston_diffusion 
    paths[,,i] <- path
    logweights[i] <- EMjointdensityMV(path,T,afun,bfun,deltat,params)-ResidualbridgejointdensityMV(path,T,heston_det,afun,bfun,deltat,params)
  }
  weights <- exp(logweights) #R can handle unnormalised weights
  mc_estimator <- mean(weights)
  if(wr) 
  {
    sampled_indices <- sample(seq(1,N,1), N, TRUE, weights)
    sampled_paths <- paths[,,sampled_indices]
    return (list(sampled_paths, log(mc_estimator)))
  }
  return(log(mc_estimator))
}

hestonSim <- function(hParams, S0,V0,T=1,deltat=1/52){
  
  mu <- hParams[1]
  kappa <- hParams[2]
  theta <- hParams[3]
  xi <- hParams[4]
  rho <- hParams[5]
  
  n <- T/deltat
  S <- rep(0,n)
  V <- rep(0,n)
  S[1] <- S0
  V[1] <- V0
  for(i in 2:n){
    BMs <- mvGauss(rho, deltat)
    Z <- BMs[1]
    W <- BMs[2]
    S[i] <- S[i-1] + mu*S[i-1]*deltat + sqrt(V[i-1])*S[i-1]*Z
    V[i] <- V[i-1] +kappa*(theta-V[i-1])*deltat + xi*sqrt(V[i-1])*W
  }
  return(list(S,V))
}

alpha_heston <- function(x,hParams){
  mu <- hParams[1]
  kappa <- hParams[2]
  theta <- hParams[3]
  return(c(mu*x[1],kappa*(theta-x[2])))
}

beta_heston <- function(x,hParams){
  xi <- hParams[4]
  rho <- hParams[5]
  var1 <- x[1]^2*x[2]
  cov <- rho*xi*x[1]*x[2]
  var2 <- xi^2*x[2]
  #print(x)
  return(matrix(c(var1,cov,cov,var2),ncol=2,byrow=TRUE))
}



loglikelihood_heston <- function(x,T=1/52,deltat=1/(2*52),N=1,params=c(0.5,0.5,0.9),afun=alpha_heston,bfun=beta_heston)
{
  if(params[5]^2>1|params[2]<=0|params[3]<=0|params[3]>1|params[4]<=0|2*params[2]*params[3]<params[4]^2){ #Assigning zero prob to impossible values of rho and alpha
    return(log(0))
  }
  logprob <- 0
  for(i in 1:(dim(x)[2]-1))
  {
    logprob <- logprob+wr_heston(x[,i],x[,i+1],T,deltat,N,params,afun,bfun)
  }
  return(logprob)
}

log_prior_heston <- function(params){
  mu <- params[1]
  kappa <- params[2]
  theta <- params[3]
  xi <- params[4]
  rho <- params[5]
  #assuming independent priors so simply add logs.
  return(dnorm(mu,log=T)+dnorm(kappa,log=T)+dnorm(theta,log=T)+dnorm(xi,log=T)+dunif(rho,min=-1,max=1,log=T))
}
mvGauss <- function(p, deltat)
{
  M <- matrix(c(deltat, p*deltat, p*deltat, deltat), nrow = 2, ncol = 2, byrow= TRUE)
  mu <- c(0,0)
  return (MASS::mvrnorm(2, mu = mu, Sigma = M))
}
EMjointdensityMV=function(x,T,afun=alphaLV_sabr,bfun=betaLV_sabr,deltat,theta)
{
  n = T/deltat
  ll=0
  for (i in 1:n){
    ll=ll+lgpdf(x[i+1,],x[i,]+afun(x[i,],theta)*deltat,bfun(x[i,],theta)*deltat)
  }
  #print(ll)
  return(ll)
}

sqrtmat=function(V)
{
  spec=svd(V)
  return(spec$u%*%diag(sqrt(spec$d))%*%t(spec$u))
}

#generates a draw from a multivariate N(m,V) distn
rmvn=function(m,V)
{
  p=length(m)
  z=rnorm(p)
  return(m+sqrtmat(V)%*%z)
}

lgpdf=function(x,m,V) #Evaluates the log of a multivariate Gaussian pdf
{
  d=length(m)
  return(-0.5*log(det(V)) -0.5*t(x-m)%*%solve(V)%*%(x-m)-0.5*d*log(2*pi))
}

