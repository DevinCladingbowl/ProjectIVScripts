library(matrixcalc)
library(emdbook)
library(MASS)
library(latex2exp)
library(Metrics)
library(MLmetrics)
library(future.apply)
library(future)
library(gmodels)

#Simulate on a finer grid than obs then thin to get obs
set.seed(1)
Xdata <- sabrSim(alpha=0.1, beta = 0.5,rho = 0.7, F0 = 100, s0=0.10, T=1, deltat=1/(10*52))
Xprice <- Xdata[[1]][1+(0:51)*10] #thin to get data points every 1/52 time units
Xvols <- Xdata[[2]][1+(0:51)*10]
par(mfrow=c(1,1))
plot(ts(Xprice,deltat=1/52,start=0))
plot(ts(Xvols,deltat=1/52,start=0))


##RUNNING SCHEME
strt <- Sys.time()
A <- sabr_latent_pmmh(Xprice,5,3,100000)
end <- Sys.time()
time <- strt-end



plot(ts(t(A[[1]]))) #Pilot run for covMat
covMat<-var(t(A[[1]])) #Then rerun

Alphs<-exp(A[[1]][1,])
Rhos <- A[[1]][2,]
plot(ts(Alphs,start=0,deltat=1),ylab=TeX("\\alpha"),xlab="Iteration",main=TeX("$\\alpha$ Trace plot"))
abline(h=0.1,col="red",lwd=3)
plot(ts(Rhos,start=0,deltat=1),ylab=TeX("\\rho"),xlab="Iteration",main=TeX("$\\rho$ Trace plot"))
abline(h=0.7,col="red",lwd=3)

hist(Alphs,freq=F,main=TeX("$\\alpha$ Histogram"),xlab=TeX("\\alpha"), breaks=100,axes=FALSE,ylim=c(0,10),mgp=c(1.5,1,0))
axis(1,at = seq(0,1,0.05),labels = TRUE,pos = 0)
axis(2,pos = 0.00)
abline(v=mean(Alphs),col="blue",lwd=4)
abline(v=0.1,col="red",lwd=4)
hist(Rhos,freq=F,main=TeX("$\\rho$ Histogram"),xlab=TeX("\\rho"), breaks=100,axes=FALSE,ylim=c(0,2),mgp=c(1.5,1,0),xlim=c(0,1))
axis(1,at = seq(0.0,1,0.05),labels = TRUE,pos = 0)
axis(2,pos = 0.0)
abline(v=mean(Rhos),col="blue",lwd=4)
abline(v=0.7,col="red",lwd=4)

midvols <- apply(A[[3]],2,mean)
lowervols <- apply(A[[3]],2,quantile,0.025)
uppervols <- apply(A[[3]],2,quantile,0.975)
plot(ts(Xvols,deltat=1/52,start=0),ylim=c(0.065,0.11),lwd=2,main="Latent path percentiles",ylab=TeX("$\\sigma_t$"))
lines(ts(lowervols,deltat=1/52,start=0),col="red",lty="dashed",lwd=2)
lines(ts(midvols,deltat=1/52,start=0),lty="dashed",lwd=2)
lines(ts(uppervols,deltat=1/52,start=0),col="blue",lty="dashed",lwd=2)

plot(acf(Alphs,lag.max=20000),main=TeX("$\\alpha$  Corellogram"),ylab=TeX("$\\rho_k"))
abline(h=0.006,col="red")
abline(h=-0.006,col="red")

#for mcmcse::multiESS()
dataf= data.frame(test[1,],test[2,],test[3,],test[4,],test[5,])


sabr_latent_pmmh <- function(Xp,N,n,M,T=1/52,deltat=1/104){ #X price series. N: #bridges. n: #intertimes.M:iterations
  
  #cov matrix estimation

  covMat <- matrix(c(0.40483540, -0.00956497, -0.00956497, 0.07031741), ncol=2,nrow=2,byrow=TRUE)
  L <- length(Xp)
  theta <- array(0,dim=c(2,M)) #top = alpha, bottom = rho
  theta[,1][1] <- log(0.1)
  theta[,1][2] <- 0.7 #initialise at ground truth
  s0 <- 0.1
  sigma_0 <- Xvols
  #plot(sigma_0)
  X<- matrix(c(Xp,sigma_0),nrow=2,byrow=TRUE)
  count <- 0
  U4D <- innovs(L,N,n,deltat) #INITIALISE INNOVATIONS ARRAY. L datapoints. N bridges per datapoint
  #can't be 1/(nL) if deltat=1/(2*52) and you're using n=3. Change last arg to deltat
  proposal_prev <- loglikelihood_latent(X,U4D,T,deltat,N,params=c(exp(theta[,1][1]),0.5, theta[,1][2])) +log_prior(theta[,1][1],theta[,1][2])
  
  ###INITIALISE LOG PROBS FOR SIGMA AND U BLOCK###
  block_logprobs <- rep(0,L)
  #block_logprobs[1] <- wr_latent_sabr(X[,1],X[,2],T,deltat,N,params=c(exp(theta[,1][1]),0.5, theta[,1][2]),U4D[1,,,],afun=alphaLV_sabr,bfun=betaLV_sabr)
    
  for (i in 2:(L-1)){
    #block_logprobs[i] <-  trunc_log(X[,(i-1):(i+1)],U4D[i-1,,,],U4D[i,,,],T,deltat,N,params=c(exp(theta[,1][1]),0.5,theta[,1][2]))
    block_logprobs[i] <-  trunc_log2(X[,(i-1):(i)],U4D[i-1,,,],T,deltat,N,params=c(exp(theta[,1][1]),0.5,theta[,1][2]))
  }
 
  block_logprobs[L] <- wr_latent_sabr(X[,(L-1)],X[,L],T,deltat,N,params=c(exp(theta[,1][1]),0.5, theta[,1][2]),U4D[L-1,,,],afun=alphaLV_sabr,bfun=betaLV_sabr)
  
  #Note only need L-1 block probs - [1] and [2] are the same in the above.
  
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
    proposal_theta <- mvrnorm(1,theta[,i-1],covMat)
    #proposal_theta <- mvn(theta[,i-1],0.1) #LOGLIKELIHOOD TAKES ENTIRE 4D ARRAY?
    target_proposal <- loglikelihood_latent(X,U4D,T,deltat,N,params=c(exp(proposal_theta[1]),0.5, proposal_theta[2])) +log_prior(proposal_theta[1],proposal_theta[2])
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
    
    
    #uncomment below for use in \sigma_0 unknown
    
    # prop_vol1 <- rnorm(1,X[,1][2],0.01) 
    # prop_U1 <- innovs(L,N,n,deltat)[1,,,] 
    # subvals1 <- X[,1]
    # subvals1[2] <- prop_vol1
    # ublockalpha1 <- wr_latent_sabr(subvals1,X[,2],T,deltat,N,params=c(exp(theta[,i][1]),0.5, theta[,i][2]),prop_U1,afun=alphaLV_sabr,bfun=betaLV_sabr)
    # if(log(runif(1)) < (ublockalpha1-block_logprobs[1])){
    #   X[,1][2] <- prop_vol1
    #   U4D[1,,,] <- prop_U1
    #   block_logprobs[1] <- ublockalpha1
    #   block_logprobs[2] <- ublockalpha1 # hacky fix given the reducndancy the block_logprobs vec
    # }
    
    prop_volL <- rnorm(1,X[,L][2],0.01) #innovation variance of 1 is huge compared to vol values! Reduce!
    prop_UL <- array(rnorm(N*(n-1)*2,0,sqrt(deltat)),dim=c(N,n-1,2))
    subvalsL <- X[,L]
    subvalsL[2] <- prop_volL
    ublockalphaL <- wr_latent_sabr(X[,L-1],subvalsL,T,deltat,N,params=c(exp(theta[,i][1]),0.5, theta[,i][2]),prop_UL,afun=alphaLV_sabr,bfun=betaLV_sabr)
    
    if(log(runif(1)) < (ublockalphaL-block_logprobs[L])){
      X[,L][2] <- prop_volL
      U4D[L-1,,,] <- prop_UL
      block_logprobs[L] <- ublockalphaL
    }
    ##THE REST##  
    for (j in 2:(L-1)){
      U3D_prev <- U4D[j-1,,,] 
      U3D_curr <- U4D[j,,,]
      prop_vol <- rnorm(1,X[,j][2],0.01) 
      if(prop_vol <0){
        next
      }
      prop_U_prev <- innovs(L,N,n,deltat)[1,,,] #PROPOSAL FOR U3D_prev
      prop_U_curr <- innovs(L,N,n,deltat)[1,,,] #PROPOSAL FOR U3D_curr
      subvals <- X[,(j-1):(j+1)]
      subvals[2,2] <- prop_vol
      proposal_logprob1 <- trunc_log2(subvals[,1:2],prop_U_prev,T,deltat,N,params=c(exp(theta[,i][1]),0.5,theta[,i][2]))
      proposal_logprob2 <- trunc_log2(subvals[,2:3],prop_U_curr,T,deltat,N,params=c(exp(theta[,i][1]),0.5,theta[,i][2]))
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
    proposal_prev <- loglikelihood_latent(X,U4D,T,deltat,N,params=c(exp(theta[,i][1]),0.5, theta[,i][2])) +log_prior(theta[,i][1],theta[,i][2])
  

  }
  print(count/(M-1))
  #return(count/(n-1))
  return(list(theta,X,vols))
}




trunc_log2 <- function(X,U0,T=1/52,deltat=1/(2*52),N=5,params=c(0.5,0.5,0.9),afun=alphaLV_sabr,bfun=betaLV_sabr)
{
  logprob <- wr_latent_sabr(X[,1],X[,2],T,deltat,N,params,U0,afun,bfun)

  return(logprob)
}

loglikelihood_latent <- function(x,U4D,T=1/52,deltat=1/(2*52),N=5,params=c(0.5,0.5,0.9),afun=alphaLV_sabr,bfun=betaLV_sabr)
{
  if(params[3]^2>1|params[1]<=0|params[3]<0){ #Assigning zero prob to impossible values of rho and alpha #rho <0 a priori
    return(log(0))
  }
  logprob <- 0
  
  for(i in 1:(dim(x)[2]-1))
  {
    U3D <- U4D[i,,,]
    logprob <- logprob+wr_latent_sabr(x[,i],x[,i+1],T,deltat,N,params,U3D,afun,bfun)
  }
  return(logprob)
}



DGbridge_latent_sabr=function(x0, xT, deltat, T, afun = alphaLV_sabr, bfun = betaLV_sabr, params=c(1,0.5,0.3),Umat) #takes sub matrix of wider 4d innovation array
{
  n = T/deltat
  d=length(x0) #d is state dimension
  x=matrix(0,ncol=d,nrow=n+1) # columns correspond to each component
  x[1,]=x0
  x[n+1,]=xT
  for(i in 1:(n-1))
  {
    # print(i)
    u <- Umat[i,]
    tm = (i-1)*deltat
    t = i*deltat
    x[i+1,] = (x[i,] + (xT-x[i,])*(deltat)/(T-tm)+sqrt((T-t)/(T-tm))*t(chol(bfun(x[i,],params)))%*%u)
    for(j in 1:d){
      if(x[i+1,j]<0.01){
        x[i+1,j]=0.01 #fudge, to avoid negative values
      }
    }
    
  }
  return(x)
}

innovs <- function(T,N,n,dt){ #T is the number of data points. (check this)
  return (array(rnorm((T-1)*N*(n-1)*2,0,sqrt(dt)),dim=c(T,N,n-1,2))) #(T-1) not T otherwise too many innovations?
}


wr_latent_sabr <- function(x0, xT, T=1, deltat=0.1, N=5, params,U3D,afun=alphaLV_sabr,bfun=betaLV_sabr,wr=FALSE) #set wr=TRUE to perform weighted resampling. Else importance weights are returned
{
  paths <- array(0,dim=c((T/deltat)+1,2,N))
  logweights <- rep(0,N)
  paths[1] <- 2
  for (i in 1:N)
  {
    Umat <- U3D[i,,]
    path <- DGbridge_latent_sabr(x0, xT, deltat, T,afun,bfun,params,Umat)
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

#Function to evaluate (log of) the joint density under proposal mechanism (modified diffusion bridge, observe all end-point components)
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


log_prior <- function(alpha,rho){
  #assuming independent priors so simply add logs.
  return(dgamma(exp(alpha),2,10,log=TRUE)+dunif(rho,min=-1,max=1,log=TRUE)+alpha) #includes log of jacobian term for prior on exp(alpha arg)
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
