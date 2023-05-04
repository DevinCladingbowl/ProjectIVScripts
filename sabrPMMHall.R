
library(MASS)
####SABR INFERENCE####
#This R script contains functions to perform inference of the SABR model using pseudo-marginal Metropolis-Hastings


###PMMH###


sabr_pmmh <- function(n, X,l){
  theta <- array(0,dim=c(2,n)) #top = alpha, bottom = rho
  theta[,1][1] <- log(0.1)
  theta[,1][2] <- 0.1
  count <- 0
  proposal_prev <- loglikelihood(X,params=c(exp(theta[,1][1]),0.5, theta[,1][2])) +log_prior(theta[,1][1],theta[,1][2])
  for (i in 2:n){
    if(i%%100==0){
      print(i)
    }
    proposal_theta <- mvn(theta[,i-1],l)
    target_proposal <- loglikelihood(X,params=c(exp(proposal_theta[1]),0.5, proposal_theta[2])) +log_prior(proposal_theta[1],proposal_theta[2])
    alpha <- target_proposal-proposal_prev
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


#####

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

alphaLV_sabr=function(x,params){
  
  return(c(0,0))
}

betaLV_sabr=function(x,params) #params= SABR alpha, beta, rho
{
  var1=x[2]*x[2]*x[1]^(2*params[2]) 
  cov=params[3]*params[1]*x[2]*x[2]*x[1]^params[2] 
  var2=(params[1]^2)*x[2]^2
  return(matrix(c(var1,cov,cov,var2),ncol=2,byrow=TRUE))
}

DGbridge_sabr=function(x0, xT, deltat, T, afun = alphaLV_sabr, bfun = betaLV_sabr, params=c(1,0.5,0.3))
{
  n = T/deltat
  d=length(x0) #d is state dimension
  x=matrix(0,ncol=d,nrow=n+1) # columns correspond to each component
  x[1,]=x0
  x[n+1,]=xT
  for(i in 1:(n-1))
  {
    tm = (i-1)*deltat
    t = i*deltat
    x[i+1,] = (x[i,] + (xT-x[i,])*(deltat)/(T-tm)+sqrt((T-t)/(T-tm))*t(chol(bfun(x[i,],params)))%*%rnorm(2,0,sqrt(deltat)))
    for(j in 1:d){
      if(x[i+1,j]<0.01){
        x[i+1,j]=0.01 #fudge, to avoid negative values
      }
    }
  }
  return(x)
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

DGbridgejointdensityMV=function(x,T,afun=alphaLV_sabr,bfun=betaLV_sabr,deltat,params)
{
  n = T/deltat
  ll=0
  x1=x[n+1,]
  for (i in 1:(n-1)){
    t=i*deltat
    ll=ll+lgpdf(x[i+1,],((T-t)*x[i,]+x1*deltat)/(T-t+deltat),((T-t)/(T-t+deltat))*deltat*bfun(x[i,],params))
  }
  return(ll)
}

wr_sabr <- function(x0, xT, T=1, deltat=0.1, N=1000, params,afun=alphaLV_sabr,bfun=betaLV_sabr,wr=FALSE) #set wr=TRUE to perform weighted resampling. Else importance weights are returned
{
  paths <- array(0,dim=c((T/deltat)+1,2,N))
  logweights <- rep(0,N)
  paths[1] <- 2
  for (i in 1:N)
  {
    path <- DGbridge_sabr(x0, xT, deltat, T,afun,bfun,params)
    paths[,,i] <- path
    logweights[i] <- EMjointdensityMV(path,T,afun,bfun,deltat,params)-DGbridgejointdensityMV(path,T,afun,bfun,deltat,params)
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

loglikelihood <- function(x,T=1/52,deltat=1/(2*52),N=1,params=c(0.5,0.5,0.9),afun=alphaLV_sabr,bfun=betaLV_sabr)
{
  if(params[3]^2>1|params[1]<=0){ #Assigning zero prob to impossible values of rho and alpha
    return(log(0))
  }
  logprob <- 0
  for(i in 1:(dim(x)[2]-1))
  {
    logprob <- logprob+wr_sabr(x[,i],x[,i+1],T,deltat,N,params,afun,bfun)
  }
  return(logprob)
}

log_prior <- function(alpha,rho){
  #assuming independent priors so simply add logs.
  return(dnorm(alpha)+log(dunif(rho,min=-1,max=1)))
}

mvn <- function(MU,l){
  M <- matrix(c(l,0,0,l),ncol=2,nrow=2,byrow=TRUE)
  return(mvrnorm(1,MU,M))
}



mvGauss <- function(p, deltat)
{
  M <- matrix(c(deltat, p*deltat, p*deltat, deltat), nrow = 2, ncol = 2, byrow= TRUE)
  mu <- c(0,0)
  return (MASS::mvrnorm(2, mu = mu, Sigma = M))
}


#Euler Maruyama for SABR. a = 0 for both. b only. Y_{n+1} = y_n + b * dW
# 0<beta<1 alpha > 0 0<rho<1
sabrSim <- function(alpha = 0.1, beta = 0.5,rho = 0.7, F0 = 100, s0=0.1, T=10, deltat=0.1) #Euler-Maruyama approximation for the SABR SDE system
{
  n <- T/deltat 
  f <- rep(0,n)
  s <- rep(0,n)
  f[1] <- F0
  s[1] <- s0
  
  for (i in 2:n)
  {
    t <- i*T/n
    BMs <- mvGauss(rho, deltat)
    Z <- BMs[1]
    W <- BMs[2]
    s[i] <- s[i-1] + alpha*s[i-1]*Z
    f[i] <- f[i-1] + s[i-1]*(f[i-1]^(beta))*W  #Z and W correlated Brownian motions. 
    #fudge for almost zero rates:
    if(s[i]<0.001){
      s[i] <- s[i-1]
    }
    if(f[i]<0.001){
      f[i] <- f[i-1]
    }
  }
  return(list(f,s))
}
