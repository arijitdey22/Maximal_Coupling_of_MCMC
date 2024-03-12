library(scales)

###=============================================================================
###Algorithm 3: (Full kernel Coupling with independent residuals)

#Target: Univariate Normal with mean = 3 and variance 5.

#rm(list = ls())

set.seed(2)

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
  
  if ( (x.cur != x.prev) & (U <=  (f.y.X / f.x.X))){
    y.cur <- x.cur
  }else{
    count.else <- count.else + 1
    accpet.y <- 0
    iter.this <- 0
    
    #while ((accpet.y == 0) & (iter.this < 500)){
    while (accpet.y == 0){
      iter.this <- iter.this + 1
      
      #proposal and acceptance of y.curl
      y.curl.cand <- rnorm(1, mean = y.cur, cand.sd)
      a.yy <- fun.a.x.y(y.cur, y.curl.cand)
      B.y <- rbinom(1, 1, min(1,a.yy))
      y.curl <- B.y * y.curl.cand + (1-B.y) * y.cur
      
      if (y.curl == y.cur){
        y.cur = y.curl
        accpet.y = 1
      }else{
        f.y.y.curl <- fun.f.x.y(y.cur, y.curl, cand.sd)
        f.x.y.curl <- fun.f.x.y(x.cur, y.curl, cand.sd)
        V <- runif(1)
        
        if (V > (f.x.y.curl / f.y.y.curl)){
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

#plotting both the chains

library(ggplot2)
library(tidyverse)
library(ggExtra)
library(gridExtra)

theme_set(theme_bw(12))

df <- data.frame(store.x, store.y) 

p <- ggplot(df, aes(x = store.x, y = store.y)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_vline( aes(xintercept = mean(store.x)), col = " darkorchid3", lwd  = 0.8, lty = "dashed") +
  geom_hline( aes(yintercept = mean(store.y)), col = " darkorchid3", lwd  = 0.8, lty = "dashed") +
  labs(x = "X", y = "Y") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))

ggMarginal(p, type = "histogram", fill = "#9bd5e8", col = "#033e52",
           binwidth = 0.5, size = 3, margins = "both")

#number of couplings
sum(store.x == store.y) / B

store.x[1:10]
store.y[1:10]   
#first two values are different only. Others are same.

