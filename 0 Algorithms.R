library(mvtnorm)
library(scales)

###=============================================================================
###Algorithm 1: (Status Quo)

#Target: Univariate Normal with mean = 3 and variance 5.

#rm(list = ls())

set.seed(100)

x.init <- 0                      #initial values
y.init <- 0

B <- 1e4                         #iterations
accept.count.x <- 0                #to track proportion of acceptance
accept.count.y <- 0

store.x <- numeric(length = B)   #defining storage
store.y <- numeric(length = B)

x.cur <- x.init                  #chain initiation
y.cur <- y.init

for (i in 1:B){
  
  #proposal
  q.xy <- rmvnorm(1, mean = c(x.cur, y.cur), sigma = 5 * diag(2))
  x.cand <- q.xy[1]
  y.cand <- q.xy[2]

  #acceptance probabilities
  a.xx <- dnorm(x.cand, 3, sqrt(5)) / dnorm(x.cur, 3, sqrt(5))
  a.yy <- dnorm(y.cand, 3, sqrt(5)) / dnorm(y.cur, 3, sqrt(5))
  
  B.x <- rbinom(1, 1, min(1, a.xx))
  B.y <- rbinom(1, 1, min(1, a.yy))
  
  #accept or reject
  x.cur <- B.x * x.cand + (1-B.x) * x.cur
  y.cur <- B.y * y.cand + (1-B.y) * y.cur
  
  accept.count.x <- accept.count.y + B.x
  accept.count.y <- accept.count.y + B.y
  
  #storing the chain
  store.x[i] <- x.cur
  store.y[i] <- y.cur
  
}

#plotting:
par(mfrow = c(1,2))
plot.ts(store.x, main = "First MCMC chain",       #plotting first chain
        xlab = "Iterations", ylab = "X")
abline(h = 3, col = "red", lwd = 2)

plot.ts(store.y, main = "Second MCMC chain",      #plotting second chain 
        xlab = "Iterations", ylab = "Y")
abline(h = 3, col = "red", lwd = 2)
par(mfrow = c(1,1))

plot(store.x, store.y, pch = 16, cex = 0.65,      #plotting both the chains
     main = "Two MCMC chains together", 
     xlab = "X", ylab = "Y")

#number of coupling
sum(store.x == store.y) / B

#acceptance prob
(acc.prob.x <- accept.count.x / B)
(acc.prob.y <- accept.count.y / B)

#parameter estimation
mean(store.x)
mean(store.y)

var(store.x)
var(store.y)


###=============================================================================
###Algorithm 2: (Maximum Coupling with independent residuals)

#Target: Univariate Normal with mean = 3 and variance 5.

#rm(list = ls())

set.seed(100)

#initial values
x.init <- runif(1)
y.init <- runif(1)

B <- 1e5                         #iterations

store.x <- numeric(length = B)   #defining storage
store.y <- numeric(length = B)

cand.sd <- 0.5                   #standard deviation of the proposal distribution
accept.count.x <- 0              #to track proportion of acceptance
accept.count.y <- 0

x.cur <- x.init                  #chain initiation
y.cur <- y.init

for (i in 1:B){
  
  #proposal of x
  x.cand <- rnorm(1, mean = x.cur, cand.sd)
  U <- runif(1)
  
  #proposal of y
  if (U <= dnorm(x.cand, y.cur, cand.sd) / dnorm(x.cand, x.cur, cand.sd)){
    y.cand <- x.cand
  }else{
    accept = 0
    while(accept == 0){
      y.curl <- rnorm(1, y.cur, cand.sd)
      V <- runif(1)
      
      if(V > dnorm(y.curl, x.cur, cand.sd) / dnorm(y.curl, y.cur, cand.sd) ){
        y.cand <- y.curl
        accept = 1
      }
    }
  }
  
  #acceptance probabilities
  a.xx <- dnorm(x.cand, 3, sqrt(5)) / dnorm(x.cur, 3, sqrt(5))
  a.yy <- dnorm(y.cand, 3, sqrt(5)) / dnorm(y.cur, 3, sqrt(5))
  
  B.x <- rbinom(1, 1, min(1, a.xx))
  B.y <- rbinom(1, 1, min(1, a.yy))
  
  #accept or reject
  x.cur <- B.x * x.cand + (1-B.x) * x.cur
  y.cur <- B.y * y.cand + (1-B.y) * y.cur
  
  accept.count.x <- accept.count.y + B.x
  accept.count.y <- accept.count.y + B.y
  
  #storing the chain
  store.x[i] <- x.cur
  store.y[i] <- y.cur
  
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

#acceptance prob
(acc.prob.x <- accept.count.x / B)
(acc.prob.y <- accept.count.y / B)

#parameter estimation
mean(store.x)
mean(store.y)

var(store.x)
var(store.y)

###----------------------------------------------------------------------------
###Recreating Figure 2:

#rm(list = ls())

set.seed(100)

#current values
x.cur <- 1/4
y.cur <- 4

N <- 1e4                         #no of replications

store.x <- numeric(length = N)   #defining storage
store.y <- numeric(length = N)

cand.sd <- 2                   #standard deviation of the proposal distribution

for (i in 1:N){
  #proposal of x
  x.cand <- rnorm(1, mean = x.cur, cand.sd)
  U <- runif(1)
  
  #proposal of y
  if (U <= dnorm(x.cand, y.cur, cand.sd) / dnorm(x.cand, x.cur, cand.sd)){
    y.cand <- x.cand
  }else{
    accept = 0
    while(accept == 0){
      y.curl <- rnorm(1, y.cur, cand.sd)
      V <- runif(1)
      
      if(V > dnorm(y.curl, x.cur, cand.sd) / dnorm(y.curl, y.cur, cand.sd) ){
        y.cand <- y.curl
        accept = 1
      }
    }
  }
  
  #acceptance probabilities
  a.xx <- dnorm(x.cand, 3, sqrt(5)) / dnorm(x.cur, 3, sqrt(5))
  a.yy <- dnorm(y.cand, 3, sqrt(5)) / dnorm(y.cur, 3, sqrt(5))
  
  B.x <- rbinom(1, 1, min(1, a.xx))
  B.y <- rbinom(1, 1, min(1, a.yy))
  
  #accept or reject
  x.acc <- B.x * x.cand + (1-B.x) * x.cur
  y.acc <- B.y * y.cand + (1-B.y) * y.cur
  
  #store
  store.x[i] <- x.acc
  store.y[i] <- y.acc
  
}

#plotting:
plot(store.x, store.y,                                
     pch = 16, cex = 0.65, col = alpha("black", 0.6),      
     main = "Coupling of MCMC from X = 0.25 and Y = 4", 
     xlab = "X", ylab = "Y")

#number of coupling
sum(store.x == store.y) / N

#creating marginals plots together

library(ggplot2)
library(tidyverse)
library(ggExtra)
library(gridExtra)

theme_set(theme_bw(12))

df <- data.frame(store.x, store.y) 

p <- ggplot(df, aes(x = store.x, y = store.y)) +
  geom_point(size = 0.8, alpha = 0.5) +
  labs(x = "X", y = "Y") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

ggMarginal(p, type = "histogram", fill = "#9bd5e8", col = "#033e52",
           binwidth = 0.5, size = 3, margins = "both")


###=============================================================================
###Algorithm 3:







