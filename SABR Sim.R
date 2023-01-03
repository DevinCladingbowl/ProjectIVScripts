#############################
###SABR FORWARD SIMULATION###
#############################

library(MASS)
mvGauss <- function(p, deltat)
{
  M <- matrix(c(deltat, p*deltat, p*deltat, deltat), nrow = 2, ncol = 2, byrow= TRUE)
  mu <- c(0,0)
  return (MASS::mvrnorm(1, mu = mu, Sigma = M))
}


# 0<beta<1 alpha > 0 0<rho<1
sabrSim <- function(alpha = 0.1, beta = 1,rho = 0.7, F0 = 100, s0=0.1, T=1, deltat=0.1) #Euler-Maruyama approximation for the SABR SDE system
{
  n <- T/deltat+1
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
    f[i] <- f[i-1] + s[i-1]*(f[i-1]^(beta))*W  #Z and W correlated Brownian motions. s[i-1] on RHS here as cond on previous!
    #fudge for almost zero rates:
    if(s[i]<0.1){
      s[i] <- s[i-1]
    }
    
  }
  return(list(f,s))
}

T <- sabrSim()
par(mfrow=c(1,2))
plot(ts(T[[1]],start=0,deltat=0.1))
plot(ts(T[[2]],start=0,deltat=0.1))

