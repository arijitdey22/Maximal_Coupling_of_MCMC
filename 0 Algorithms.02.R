library(scales)

###=============================================================================
###Algorithm 3: (Full kernel Coupling with independent residuals)

#Target: Univariate Normal with mean = 3 and variance 5.

rm(list = ls())

set.seed(100)

B <- 1e4                         #iterations

store.x <- numeric(length = B)   #defining storage
store.y <- numeric(length = B)

accept.count.x <- 0              #to track proportion of acceptance

#-- --
#necessary functions

fun.a.x.y <- function(x.1,y.1){
  min(1, dnorm(y.1, 3, sqrt(5)) / dnorm(x.1, 3, sqrt(5)))
}

fun.q.x.y <- function(x.2, y.2, cand.sd){
  dnorm(y.2, x.2, cand.sd)
}

fun.f.x.y <- function(x.3, y.3, cand.sd){
  fun.a.x.y(x.3,y.3) * fun.q.x.y(x.3, y.3, cand.sd)
}

#-- --

cand.sd <- 1                   #standard deviation of the proposal distribution

#initial values
x.cur <- runif(1)
y.cur <- runif(1)

count.else <- 0

for (i in 1:B){
  
  #proposal of x
  x.cand <- rnorm(1, mean = x.cur, cand.sd) 
  
  #accepting x
  a.xx <- fun.a.x.y(x.cur, x.cand)
  B.x <- rbinom(1, 1, min(1, a.xx))
  x.prev <- x.cur
  x.cur <- B.x * x.cand + (1-B.x) * x.cur
  
  accept.count.x <- accept.count.x + B.x
  
  #accepting y
  f.x.X <- fun.f.x.y(x.prev, x.cur, cand.sd)
  f.y.X <- fun.f.x.y(y.cur, x.cur, cand.sd)
  U <- runif(1)
  
  if ( (x.cur != x.prev) & (U <=  (f.x.X / f.y.X))){
    y.cur <- x.cur
  }else{
    count.else <- count.else + 1
    accpet.y <- 0
    iter.this <- 0
    
    #while ((accpet.y == 0) & (iter.this < 500)){
    while (accpet.y == 0){
      iter.this <- iter.this + 1
      
      y.curl <- rnorm(1, mean = y.cur, cand.sd)
      
      if (y.curl == y.cur){
        y.cur = y.curl
        accpet.y = 1
      }else{
        f.y.y.curl <- fun.f.x.y(y.cur, y.curl, cand.sd)
        f.x.y.curl <- fun.f.x.y(x.cur, y.curl, cand.sd)
        V <- runif(1)
        
        if (V < (f.x.y.curl / f.y.y.curl)){
          y.cur = y.curl
          accpet.y = 1
        }
      }
      
    }
    
  }
  
  #storing the chain
  store.x[i] <- x.cur
  store.y[i] <- y.cur
  
  #aesthetics
  if(i %% 100 == 0){
    print(paste0(i, " iterations completed."))
  }
  
}

#plotting:
par(mfrow = c(1,2))
plot.ts(store.x, main = "First MCMC chain",                                #plotting first chain
        xlab = "Iterations", ylab = "X")
abline(h = 3, col = "red", lwd = 2)

plot.ts(store.y, main = "Second MCMC chain",                               #plotting second chain 
        xlab = "Iterations", ylab = "Y")
abline(h = 3, col = "red", lwd = 2)
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
mean(store.x)
mean(store.y)

var(store.x)
var(store.y)

