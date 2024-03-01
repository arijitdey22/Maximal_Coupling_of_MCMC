library(scales)
library(mvtnorm)

###=============================================================================
###Algorithm 3: (Full kernel Coupling with independent residuals)

#Target: Multivariate Normal with mean = rep(2,d) and sigma:
    #Case I: 2 * diag(d)
    #Case II: CAR structure. We do this multiple times for varying rho.

rm(list = ls())

set.seed(100)

d <- 2                           #dimension 
B <- 1e4                         #iterations
target.sigma <- 2 * diag(d)

# foo <- matrix(rep(1:d, each = d), ncol = d)
# target.sigma <- 0.5 ^ abs(foo - t(foo))

store.x <- matrix(0, B, d)       #defining storage
store.y <- matrix(0, B, d)       

accept.count.x <- 0              #to track proportion of acceptance

#-- --
#necessary functions

fun.a.x.y <- function(x.1,y.1, d, target.sigma){
  min(1, dmvnorm(y.1, mean = rep(2,d), sigma = target.sigma) / dmvnorm(x.1, mean = rep(2,d), sigma = target.sigma))
}

fun.q.x.y <- function(x.2, y.2, cand.sd){
  dmvnorm(y.2, x.2, cand.sd^2 * diag(d))
}

fun.f.x.y <- function(x.3, y.3, cand.sd, d, target.sigma){
  fun.a.x.y(x.3,y.3, d, target.sigma) * fun.q.x.y(x.3, y.3, cand.sd)
}

#-- --

cand.sd <- 1                   #standard deviation of the proposal distribution

#initial values
x.cur <- runif(d)
y.cur <- runif(d)

count.else <- 0

for (i in 1:B){
  
  #proposal of x
  x.cand <- rmvnorm(1, mean = x.cur, sigma = cand.sd^2 * diag(d))
  
  #accepting x
  a.xx <- fun.a.x.y(x.cur, x.cand, d, target.sigma)
  B.x <- rbinom(1, 1, min(1, a.xx))
  x.prev <- x.cur
  x.cur <- B.x * x.cand + (1-B.x) * x.cur
  
  accept.count.x <- accept.count.x + B.x
  
  #accepting y
  f.x.X <- fun.f.x.y(x.prev, x.cur, cand.sd, d, target.sigma)
  f.y.X <- fun.f.x.y(y.cur, x.cur, cand.sd, d, target.sigma)
  U <- runif(1)
  
  if ( (sum(x.cur == x.prev) == 0) & (U <=  (f.x.X / f.y.X))){
    y.cur <- x.cur
  }else{
    count.else <- count.else + 1
    accpet.y <- 0
    iter.this <- 0
    
    #while ((accpet.y == 0) & (iter.this < 500)){
    while (accpet.y == 0){
      iter.this <- iter.this + 1
      
      y.curl <- rmvnorm(1, mean = y.cur, sigma = cand.sd^2 * diag(d))
      
      if (sum(y.curl != y.cur) == 0){
        y.cur = y.curl
        accpet.y = 1
      }else{
        f.y.y.curl <- fun.f.x.y(y.cur, y.curl, cand.sd, d, target.sigma)
        f.x.y.curl <- fun.f.x.y(x.cur, y.curl, cand.sd, d, target.sigma)
        V <- runif(1)
        
        if (V < (f.x.y.curl / f.y.y.curl)){
          y.cur = y.curl
          accpet.y = 1
        }
      }
      
    }
    
  }
  
  #storing the chain
  store.x[i, ] <- x.cur
  store.y[i, ] <- y.cur
  
  #aesthetics
  if(i %% 100 == 0){
    print(paste0(i, " iterations completed."))
  }
  
}

#plotting:
par(mfrow = c(1,2))
plot.ts(store.x[,1], main = "First MCMC chain",                                #plotting first chain
        xlab = "Iterations", ylab = "X")
abline(h = 2, col = "red", lwd = 2)

plot.ts(store.y[,1], main = "Second MCMC chain",                               #plotting second chain 
        xlab = "Iterations", ylab = "Y")
abline(h = 2, col = "red", lwd = 2)
par(mfrow = c(1,1))

plot(store.x[store.x != store.y], store.y[store.x != store.y],             #plotting both the chains
     pch = 16, cex = 0.65, col = alpha("black", 0.6),      
     main = "Two MCMC chains together", 
     xlab = "X", ylab = "Y")
points(store.x[store.x == store.y], store.y[store.x == store.y],
       col = alpha("red", 0.8), pch = 16, cex = 0.65)

#number of coupling
sum(store.x == store.y) / B

#acceptance prob of x
(acc.prob.x <- accept.count.x / B)

#number of elses
count.else / B

#parameter estimation
colMeans(store.x)
colMeans(store.y)

apply(store.x, 2, var)
apply(store.y, 2, var)

