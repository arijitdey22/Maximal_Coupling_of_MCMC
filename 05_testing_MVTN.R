library(scales)
library(mvtnorm)

###=============================================================================
###Testing the time taken for Algo 3

#Target: Multivariate Normal with mean = rep(2,d) and 
#                                 sigma: CAR structure.

rm(list = ls())

##--------------------------------------
#necessary functions for algorithms

my.dmvnorm <- function(X, mu, sigma.inv, cons){
  
  deviation <- X - mu
  exp.term <- as.numeric(t(deviation) %*% sigma.inv %*% deviation)
  density <- cons * exp( -0.5 * exp.term )
  return(density)
  
}

fun.a.x.y <- function(x.1,y.1, mu, sigma.inv, cons){
  min(1, my.dmvnorm(y.1, mu, sigma.inv, cons) / my.dmvnorm(x.1, mu, sigma.inv, cons))
}

fun.q.x.y <- function(x.2, y.2, cand.sd){
  prod(dnorm(y.2, x.2, cand.sd))
}

fun.f.x.y <- function(x.3, y.3, cand.sd, mu, sigma.inv, cons){
  fun.a.x.y(x.3,y.3, mu, sigma.inv, cons) * fun.q.x.y(x.3, y.3, cand.sd)
}

##--------------------------------------

d <- 30
cand.sd <- 1
y.accept.alg.3 <- 0

target.mean <- rep(2,d)                               #target distribution

foo <- matrix(rep(1:d, each = d), ncol = d)         #case 2
target.sigma <- 0.5 ^ abs(foo - t(foo))

  target.sigma.inv <- solve(target.sigma)
  cons <- 1 / sqrt( det(2*pi *target.sigma) )

x.cur.3 <- runif(d)                            #initial values
y.cur.3 <- runif(d)

y.accept.alg.3 <- 0                            #marker variable for acceptance
count.iter.3 <- 0                              #count of the first meeting
while.loop.count <- 0                          #count of the while loop

while (y.accept.alg.3 == 0){
  
  count.iter.3 <- count.iter.3 + 1
  
  #proposal of x
  x.cand.3 <- x.cur.3 + cand.sd *  rnorm(d)
  
  #accepting x
  a.xx.3 <- fun.a.x.y(x.cur.3, x.cand.3, target.mean, target.sigma.inv, cons)
  B.x.3 <- rbinom(1, 1, min(1, a.xx.3))
  x.prev.3 <- x.cur.3
  x.cur.3 <- B.x.3 * x.cand.3 + (1-B.x.3) * x.cur.3
  
  #accepting y
  f.x.X.3 <- fun.f.x.y(x.prev.3, x.cur.3, cand.sd, target.mean, target.sigma.inv, cons)
  f.y.X.3 <- fun.f.x.y(y.cur.3, x.cur.3, cand.sd, target.mean, target.sigma.inv, cons)
  U <- runif(1)
  
  if ( (sum(x.cur.3 == x.prev.3) == 0) & (U <=  (f.x.X.3 / f.y.X.3))){
    y.cur.3 <- x.cur.3
    y.accept.alg.3 <- 1
  }else{
    
    accpet.y.3 <- 0
    
    while (accpet.y.3 == 0){
      
      while.loop.count <- while.loop.count + 1
      
      y.curl.3.cand <- y.cur.3 + cand.sd * rnorm(d)
      a.yy.3 <- fun.a.x.y(y.cur.3, y.curl.3.cand, target.mean, target.sigma.inv, cons)
      B.y.3 <- rbinom(1, 1, min(1, a.yy.3))
      y.curl.3 <- B.y.3 * y.curl.3.cand + (1-B.y.3) * y.cur.3
      
      if (sum(y.curl.3 != y.cur.3) == 0){
        y.cur.3 = y.curl.3
        accpet.y.3 = 1
      }else{
        f.y.y.curl.3 <- fun.f.x.y(y.cur.3, y.curl.3, cand.sd, target.mean, target.sigma.inv, cons)
        f.x.y.curl.3 <- fun.f.x.y(x.cur.3, y.curl.3, cand.sd, target.mean, target.sigma.inv, cons)
        V <- runif(1)
        
        if (V < (f.x.y.curl.3 / f.y.y.curl.3)){
          y.cur.3 = y.curl.3
          accpet.y.3 = 1
          
        }
      }
      
    }
    
  }
  
}

while.loop.count
