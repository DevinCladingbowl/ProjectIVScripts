library(matrixcalc)
library(emdbook)
library(MASS)
library(latex2exp)
library(Metrics)
library(MLmetrics)
library(future.apply)
library(future)
library(gmodels)



set.seed(3005)
parameters <- c(1,5,0.15,0.1,0.7) #mu kappa theta xi rho
Hdata <- hestonSim(hParams=parameters,S0=100,V0=0.1,T=1,deltat=1/(10*52))
Hprice <- Hdata[[1]][1+(0:51)*10] #thin to get data points every 1/52 time units
Hvols <- Hdata[[2]][1+(0:51)*10]

plot(ts(Hprice,deltat=1/52,start=0))
plot(ts(Hvols,deltat=1/52,start=0))
strt <- Sys.time()
H <- heston_latent_pmmh(Hprice,5,3,100000)
end <- Sys.time()
time <- strt-end
print(time)

mu<-H[[1]][1,]
kappa <-exp(H[[1]][2,])
theta <- exp(H[[1]][3,])
xi <- exp(H[[1]][4,])
rho <- H[[1]][5,]

plot(ts(mu,start=0,deltat=1),ylab=TeX("\\mu"),xlab="Iteration",main=TeX("$\\mu$ Trace plot"))
plot(ts(kappa,start=0,deltat=1),ylab=TeX("\\kappa"),xlab="Iteration",main=TeX("$\\kappa$ Trace plot"))
plot(ts(theta,start=0,deltat=1),ylab=TeX("\\theta"),xlab="Iteration",main=TeX("$\\theta$ Trace plot"))
plot(ts(xi,start=0,deltat=1),ylab=TeX("\\xi"),xlab="Iteration",main=TeX("$\\xi$ Trace plot"))
plot(ts(rho,start=0,deltat=1),ylab=TeX("\\rho"),xlab="Iteration",main=TeX("$\\rho$ Trace plot"))


hist(mu,freq=F,main=TeX("$\\mu$ Histogram"),xlab=TeX("\\mu"), axes=FALSE,breaks=100,ylim=c(0,3),mgp=c(1.5,1,0),xlim=c(0.5,2))
axis(1,at = seq(0.5,2,0.5),labels = TRUE,pos = 0)
axis(2,pos = 0.5)
mean(test[1,])
abline(v=mean(mu),col="blue",lwd=4)
abline(v=parameters[1],col="red",lwd=4)

hist(kappa,freq=F,main=TeX("$\\kappa$ Histogram"),xlab=TeX("\\kappa"), axes=FALSE,breaks=100,ylim=c(0,3),mgp=c(1.5,1,0),xlim=c(0,6))
exp(mean(test[2,]))
axis(1,at = seq(0,6,0.5),labels = TRUE,pos = 0)
axis(2,pos = 0)
abline(v=mean(kappa),col="blue",lwd=4)
abline(v=parameters[2],col="red",lwd=4)

hist(theta,freq=F,main=TeX("$\\theta$ Histogram"),xlab=TeX("\\theta"),axes=FALSE, breaks=100,ylim=c(0,5),mgp=c(1.5,1,0),xlim=c(0.05,0.7))
axis(1,at = seq(0.05,0.7,0.05),labels = TRUE,pos = 0)
axis(2,pos = 0.05)
exp(mean(test[3,]))
abline(v=mean(theta),col="blue",lwd=4)
abline(v=parameters[3],col="red",lwd=4)

hist(xi,freq=F,main=TeX("$\\xi$ Histogram"),xlab=TeX("\\xi"),axes=FALSE, breaks=100,ylim=c(0,5),mgp=c(1.5,1,0),xlim=c(0,0.45))
axis(1,at = seq(0,0.45,0.05),labels = TRUE,pos = 0)
axis(2,pos = 0)
exp(mean(test[4,]))
abline(v=mean(xi),col="blue",lwd=4)
abline(v=parameters[4],col="red",lwd=4)

hist(rho,freq=F,main=TeX("$\\rho$ Histogram"),xlab=TeX("\\rho"), breaks=100,ylim=c(0,2),mgp=c(1.5,1,0),xlim=c(0,1),axes=FALSE)
axis(1,at = seq(0,1,0.5),labels = TRUE,pos = 0)
axis(2,pos = 0)
mean(test[5,])
abline(v=mean(rho),col="blue",lwd=4)
abline(v=parameters[5],col="red",lwd=4)


lowervols <- apply(H[[3]],MARGIN=2,FUN=function(V){return(quantile(V,probs=0.025))})
uppervols <- apply(H[[3]],MARGIN=2,FUN=function(V){return(quantile(V,probs=0.975))})
midvols <- apply(H[[3]],MARGIN=2,FUN=function(V){return(quantile(V,probs=0.5))})
plot(ts(Hvols,deltat=1/52,start=0),ylim=c(0,0.3),lwd=2,main="Latent path percentiles",ylab=TeX("$\\nu_t$"))
lines(ts(lowervols,deltat=1/52,start=0),col="red",lty="dashed",lwd=2)
lines(ts(midvols,deltat=1/52,start=0),lty="dashed",lwd=2)
lines(ts(uppervols,deltat=1/52,start=0),col="blue",lty="dashed",lwd=2)

dataf= data.frame(mu,kappa,theta,xi,rho)

plot(acf(xi,lag.max=5000),main=TeX("$\\xi$  Corellogram"),ylab=TeX("$\\rho_k"))
abline(h=0.006,col="red")
abline(h=-0.006,col="red")
dataf= data.frame(test[1,],test[2,],test[3,],test[4,],test[5,])


heston_latent_pmmh <- function(Xp,N,n,M,T=1/52,deltat=1/104){ #X price series. N: #bridges. n: #intertimes.M:iterations

  L <- length(Xp)
  theta <- array(0,dim=c(5,M)) #top = alpha, bottom = rho
  
  theta[,1][1] <- 0.1
  theta[,1][2] <- log(0.1)
  theta[,1][3] <- log(0.1)
  theta[,1][4] <- log(0.1)
  theta[,1][5] <- 0.7
  s0 <- 0.1
  sigma_0 <- Hvols
  #plot(sigma_0)
  X<- matrix(c(Xp,sigma_0),nrow=2,byrow=TRUE)
  count <- 0
  U4D <- innovs(L,N,n,deltat) #INITIALISE INNOVATIONS ARRAY. L datapoints. N bridges per datapoint

  input_params <- vector()
  input_params[1] <- theta[,1][1]
  input_params[2:4] <- exp(theta[,1][2:4])
  input_params[5] <- theta[,1][5]
  proposal_prev <- loglikelihood_latent_heston(X,U4D,T,deltat,N,params=input_params) +log_prior_heston(theta[,1])
  
  ###INITIALISE LOG PROBS FOR SIGMA AND U BLOCK###
  block_logprobs <- rep(0,L)
  block_logprobs[1] <- wr_latent_heston(X[,1],X[,2],T,deltat,N,params=input_params,U4D[1,,,],afun=alpha_heston,bfun=beta_heston)
  
  for (i in 2:(L-1)){
    #block_logprobs[i] <-  trunc_log(X[,(i-1):(i+1)],U4D[i-1,,,],U4D[i,,,],T,deltat,N,params=c(exp(theta[,1][1]),0.5,theta[,1][2]))
    block_logprobs[i] <-  trunc_log2_heston(X[,(i-1):(i)],U4D[i-1,,,],T,deltat,N,params=input_params)
  }
  
  block_logprobs[L] <- wr_latent_heston(X[,(L-1)],X[,L],T,deltat,N,params=input_params,U4D[L-1,,,],afun=alpha_heston,bfun=beta_heston)
  
  #ARRAY TO STORE ALL VOLATILITY PATHS
  vols <- array(dim=c(M,length(sigma_0)))
  vols[1,] <- sigma_0
  
  ##PMMH ITERATIONS##
  for (i in 2:M){
    if(i%%100==0){
      #plot(ts(X[2,],start=0,deltat=1/52))
      print(i)
      
    }
    
    ##########PARAMETER UPDATE#############
    #######################################
    proposal_theta <- rnorm(5,theta[,i-1],0.1)
    input_proposal <- vector()
    input_proposal[1] <- proposal_theta[1]
    input_proposal[2:4] <- exp(proposal_theta[2:4])
    input_proposal[5] <- proposal_theta[5]
    
    target_proposal <- loglikelihood_latent_heston(X,U4D,T,deltat,N,params=input_proposal) +log_prior_heston(proposal_theta)
    alpha1 <- target_proposal-proposal_prev
    if (log(runif(1))<alpha1){
      theta[,i] <- proposal_theta
      proposal_prev <- target_proposal
      count <- count + 1
    }
    else {
      theta[,i] <- theta[,i-1]
    }
    
    ############################################## 
    #######VOLATILITY AND U3D BLOCK UPDATE########
    ##ENDPOINTS#  
    input_params <- vector()
    input_params[1] <- theta[,i][1]
    input_params[2:4] <- exp(theta[,i][2:4])
    input_params[5] <- theta[,i][5]
    
    prop_volL <- rnorm(1,X[,L][2],0.01) #innovation variance of 1 is huge compared to vol values! Reduce!
    prop_UL <- array(rnorm(N*(n-1)*2,0,sqrt(deltat)),dim=c(N,n-1,2))
    subvalsL <- X[,L]
    subvalsL[2] <- prop_volL
    ublockalphaL <- wr_latent_heston(X[,L-1],subvalsL,T,deltat,N,params=input_params,prop_UL,afun=alpha_heston,bfun=beta_heston)
    
    if(log(runif(1)) < (ublockalphaL-block_logprobs[L])){
      X[,L][2] <- prop_volL
      U4D[L-1,,,] <- prop_UL
      block_logprobs[L] <- ublockalphaL
    }
    ##THE REST##  
    for (j in 2:(L-1)){
      U3D_prev <- U4D[j-1,,,] 
      U3D_curr <- U4D[j,,,]
      #prop_vol <- runif(1,0,1) #PLACEHOLDER PROPOSAL OF VOL
      prop_vol <- rnorm(1,X[,j][2],0.01)
      while(prop_vol<0.01){prop_vol <- rnorm(1,X[,j][2],0.01)} #stop chain getting sticky with low vol values
      prop_U_prev <- innovs(L,N,n,deltat)[1,,,] #PROPOSAL FOR U3D_prev
      prop_U_curr <- innovs(L,N,n,deltat)[1,,,] #PROPOSAL FOR U3D_curr
      subvals <- X[,(j-1):(j+1)]
      subvals[2,2] <- prop_vol
      proposal_logprob1 <- trunc_log2_heston(subvals[,1:2],prop_U_prev,T,deltat,N,params=input_params)
      proposal_logprob2 <- trunc_log2_heston(subvals[,2:3],prop_U_curr,T,deltat,N,params=input_params)
      alpha2 <- proposal_logprob1 + proposal_logprob2 - block_logprobs[j] - block_logprobs[j+1]
      if(log(runif(1))<alpha2){
        X[,j][2] <- prop_vol
        U4D[j-1,,,] <- prop_U_prev
        U4D[j,,,] <- prop_U_curr
        block_logprobs[j] <- proposal_logprob1
        block_logprobs[j+1] <- proposal_logprob2
      }
      else{
        #print("REJECTION!")
      }
    }
    
    vols[i,] <- X[2,]
    proposal_prev <- loglikelihood_latent_heston(X,U4D,T,deltat,N,params=input_params) +log_prior_heston(theta[,i])
  }
  print(count/(M-1))
  #return(count/(n-1))
  return(list(theta,X,vols))
}


trunc_log2_heston <- function(X,U0,T=1/52,deltat=1/(2*52),N=5,params=c(0.5,0.5,0.9),afun=alpha_heston,bfun=beta_heston)
{
  logprob_prev <- wr_latent_heston(X[,1],X[,2],T,deltat,N,params,U0,afun,bfun)
  return(logprob_prev)
}

wr_latent_heston <- function(x0, xT, T=1, deltat=0.1, N=5, params,U3D,afun=alpha_heston,bfun=beta_heston,wr=FALSE) #set wr=TRUE to perform weighted resampling. Else importance weights are returned
{
  paths <- array(0,dim=c((T/deltat)+1,2,N))
  logweights <- rep(0,N)
  paths[1] <- 2
  #print(xT)
  for (i in 1:N)
  {
    Umat <- U3D[i,,]
    heston_det <- heston_det_func(params,x0,T,deltat)
    path <- heston_residual_latent(x0, xT, heston_det,deltat, T,afun,bfun,params,Umat)
    paths[,,i] <- path
    logweights[i] <- EMjointdensityMV(path,T,afun,bfun,deltat,params)-ResidualbridgejointdensityMV(path,T,heston_det,afun,bfun,deltat,params)
  }
  max_logweights <- max(logweights)
  weights <- exp(logweights - max_logweights) # subtract the maximum value from the log weights
  mc_estimator <- sum(weights) / N
  mc_estimator <- log(mc_estimator) + max_logweights
  if(wr) 
  {
    sampled_indices <- sample(seq(1,N,1), N, TRUE, weights)
    sampled_paths <- paths[,,sampled_indices]
    return (list(sampled_paths, mc_estimator))
  }
  #print(mc_estimator)
  return(mc_estimator)
}


heston_residual_latent <- function(x0,xT,eta,deltat,T,afun=alpha_heston,bfun=beta_heston,params,Umat){ #bfun is just the heston mdb #afun is bit different
  n = T/deltat
  d=length(x0) #d is state dimension
  x=matrix(0,ncol=d,nrow=n+1) # columns correspond to each component
  rb <- matrix(0,ncol=d,nrow=n+1)
  x[1,]=x0
  x[n+1,]=xT
  nT <- eta[n+1,]
  for(i in 1:(n-1))
  {
    u <- Umat[i,]
    delta_eta <- (eta[i+1,]-eta[i,])/deltat   #The eta defined in thesis and 'Improved bridge constructs'
    tm = (i-1)*deltat
    t = i*deltat
    #print(x[i,])
    x[i+1,] = (x[i,]+ delta_eta*deltat +(xT-x[i,]-nT+eta[i,])*(deltat)/(T-tm)+sqrt((T-t)/(T-tm))*t(chol(bfun(x[i,],params)))%*%u)
    for(j in 1:d){
      if(x[i+1,j]<0.01){
        x[i+1,j]=0.01 #fudge, to avoid negative values
      }
    }
  }
  return(x)
}

loglikelihood_latent_heston <- function(x,U4D,T=1/52,deltat=1/(2*52),N=5,params,afun=alpha_heston,bfun=beta_heston)
{
  if(params[5]^2>1|params[2]<=0|params[3]<=0|params[3]>1|params[1]==0|params[5]<=0|params[4]<=0|2*params[2]*params[3]<params[4]^2){ #Assigning zero prob to impossible values of rho and alpha
    return(log(0))
  }
  logprob <- 0
  for(i in 1:(dim(x)[2]-1))
  {
    U3D <- U4D[i,,,]
    logprob <- logprob+wr_latent_heston(x[,i],x[,i+1],T,deltat,N,params,U3D,afun,bfun)
    #print(i)
  }
  return(logprob)
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

log_prior_heston <- function(params){
  mu <- params[1]
  kappa <- params[2]
  theta <- params[3]
  xi <- params[4]
  rho <- params[5]
  #assuming independent priors so simply add logs.
  return(dnorm(mu,log=T)+dnorm(kappa,log=T)+dnorm(theta,log=T)+dnorm(xi,log=T)+dunif(rho,min=-1,max=1,log=T))
}

mvn <- function(MU,l){
  M <- matrix(c(l,0,0,l),ncol=2,nrow=2,byrow=TRUE)
  return(mvrnorm(1,MU,M))
}

mvGauss <- function(p, deltat)
{
  M <- matrix(c(deltat, p*deltat, p*deltat, deltat), nrow = 2, ncol = 2, byrow= TRUE)
  mu <- c(0,0)
  return (MASS::mvrnorm(1, mu = mu, Sigma = M))
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

innovs <- function(T,N,n,dt){ #T is the number of data points. (check this)
  return (array(rnorm((T-1)*N*(n-1)*2,0,sqrt(dt)),dim=c(T,N,n-1,2))) #(T-1) not T otherwise too many innovations?
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

em <- beta_heston(c(123,0.15),c(0.1,0.1,0.1,0.1,0.990))
lgpdf(c(101,0.11),(c(100,0.10)+alpha_heston(c(123,0.15),c(0.1,0.1,0.1,0.1,0.7)))*1/104,em*1/104)

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

