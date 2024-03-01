# #Target: Normal distribution with parameters 2 and 5
# 
# #-- -- -- 
# #necessary functions
# 
# T.x.y.dens <- function(x, y, cand.sd){
#   dnorm(y, x, cand.sd)
# }
# 
# T.x..samp <- function(x, cand.sd, n = 1){
#   rnorm(n, x, cand.sd)
# }
# 
# w.x.y <- function(x, y, cand.sd){
#   lamb.x.y <- 2 / ( T.x.y.dens(x,y,cand.sd) + T.x.y.dens(y,x,cand.sd) )
#   dnorm(x,2,sqrt(5)) * T.x.y.dens(x,y,cand.sd) * lamb.x.y
# }
# 
# B <- 5e4
# k <- 10
# 
# cand.sd <- 1
# store.x <- numeric(length = B)
# 
# x.cur <- rnorm(1)
# 
# for (i in 1:B){
#   
#   set.seed(i)
#   
#   y.all <- T.x..samp(x.cur,cand.sd,k)
#   w.y.all.x <- w.x.y(y.all, x.cur, cand.sd)
#   
#   y <- sample(y.all, 1, prob = w.y.all.x)
#   
#   x.star.all <- T.x..samp(y,cand.sd,k-1)
#   x.star.all <- c(x.star.all, x.cur)
#   
#   g.hr <- min(1,sum(w.x.y(y.all,x.cur,cand.sd)) / sum(w.x.y(x.star.all,y,cand.sd)))
#   
#   B <- rbinom(1,1,g.hr)
#   x.cur <- y * B + x.cur * (1-B)
#   
#   store.x[i] <- x.cur
# }
# 
# #-- -- -- 
# #parameter estimation
# mean(store.x)
# var(store.x)
# 
# #plotting
# plot(store.x, type = "l",
#      xlab = "Iterations", ylab = "X",
#      main = paste0("MTM MCMC; Mean(X) = ", round(mean(store.x), 3)))
# abline(a = mean(store.x), b = 0, lwd = 2, col = "red")

#=============================================================

#Target: Bi-variate normal distribution with parameters
# c(1,2) and covariance matrix matrix(c(5,-2,-2,3), ncol = 2)

rm(list = ls())

library(mvtnorm)

#-- -- -- 
#necessary functions

mu.tar <- c(1,2)
sigma.tar <- matrix(c(5,-2,-2,3), ncol = 2)

T.x.y.dens <- function(x, y, cand.sd){
  dmvnorm(y, x, cand.sd^2 * diag(2))
}

T.x..samp <- function(x, cand.sd, n = 1){
  rmvnorm(n, x, cand.sd^2*diag(2))
}

w.x.y <- function(x, y, cand.sd){
  lamb.x.y <- 2 / ( T.x.y.dens(x,y,cand.sd) + T.x.y.dens(y,x,cand.sd) )
  dmvnorm(x,mu.tar,sigma.tar) * T.x.y.dens(x,y,cand.sd) * lamb.x.y
}

B <- 1e4
k <- 10

cand.sd <- 1
store.x <- matrix(0, ncol = 2, nrow = B)

x.cur <- rnorm(2)

for (i in 1:B){
  
  set.seed(i)
  
  y.all <- T.x..samp(x.cur,cand.sd,k)
  w.y.all.x <- apply(y.all, 1, w.x.y, y = x.cur, cand.sd = cand.sd)

  y <- y.all[sample(1:nrow(y.all), 1, prob = w.y.all.x), ]
  
  x.star.all <- T.x..samp(y,cand.sd,k-1)
  x.star.all <- rbind(x.star.all, x.cur)
  
  num <- sum( apply(y.all, 1, w.x.y, y = x.cur, cand.sd = cand.sd) )
  den <- sum(apply(x.star.all, 1, w.x.y, y = y, cand.sd = cand.sd))
  
  g.hr <- min(1, num  / den )
  
  B <- rbinom(1,1,g.hr)
  x.cur <- y * B + x.cur * (1-B)
  
  store.x[i,] <- x.cur
  print(i)
}

#-- -- -- 
#parameter estimation
colMeans(store.x)
cov(store.x)

#plotting
par(mfrow = c(1,2))

plot(store.x[,1], type = "l",
     xlab = "Iterations for", ylab = "X",
     main = paste0("MTM MCMC; Mean(X) = ", round(mean(store.x[,1]), 3)),
     sub = "the first component of the variable")
abline(a = mean(store.x[,1]), b = 0, lwd = 2, col = "red")

plot(store.x[,2], type = "l",
     xlab = "Iterations for", ylab = "Y",
     main = paste0("MTM MCMC; Mean(Y) = ", round(mean(store.x[,2]), 3)),
     sub = "the second component of the variable")
abline(a = mean(store.x[,2]), b = 0, lwd = 2, col = "red")

par(mfrow = c(1,1))
