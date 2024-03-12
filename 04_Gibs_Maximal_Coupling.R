# ###============================================================================
# ### Parameter estimate of N(mu, sigma) form data in Bayesian set-up
# 
# rm(list = ls())
# 
# library(invgamma)
# set.seed(1)
# 
# #-- -- 
# #data generations
# 
# mu <- 5
# sigmaSq <- 1
# 
# Y <- rnorm(500, mean = mu, sd = sqrt(sigmaSq))
# 
# #-- -- 
# #full conditional distributions
# #of mu given sigma
# 
# mu.update.samp <- function(Y, sigmaSq, mu0, sigmaSq0){
#   mu.posvar <- 1 / (length(Y) / sigmaSq + 1 / sigmaSq0)
#   mu.posmean <- (sum(Y) / sigmaSq + mu0 / sigmaSq0) * mu.posvar
#   out <- mu.posmean + sqrt(mu.posvar) * rnorm(1)
#   return(out)
# }
# 
# mu.update.den <- function(mu.fn, Y, sigmaSq, mu0, sigmaSq0){
#   mu.posvar <- 1 / (length(Y) / sigmaSq + 1 / sigmaSq0)
#   mu.posmean <- (sum(Y) / sigmaSq + mu0 / sigmaSq0) * mu.posvar
#   out <- dnorm(mu.fn, mu.posmean, sqrt(mu.posvar))
#   return(out)
# }
# 
# #of sigma given mu
# 
# sigmaSq.update.samp <- function(Y, mu, a, b){
#   sigmaSq.inv <- rgamma(1, shape = a + length(Y) / 2,
#                         rate = b + sum((Y- mu)^2) / 2)
#   out <- 1 / sigmaSq.inv
#   return(out)
# }
# 
# sigmaSq.update.den <- function(sigmaSq.fn, Y, mu, a, b){
#   if (sigmaSq.fn <= 0){
#     out <- 0
#   }else{
#     out <- dinvgamma(sigmaSq.fn, shape = a + length(Y) / 2,
#                           rate = b + sum((Y- mu)^2) / 2)
#   }
# 
#   return(out)
# }
# 
# #of the p kernel
# 
# p.samp <- function(mu.g, sigmasq.g, Y, mu0, sigmaSq0, a, b){
#   mu.n <- mu.update.samp(Y, sigmasq.g, mu0, sigmaSq0)
#   sigmasq.n <- sigmaSq.update.samp(Y, mu.n, a, b)
#   return(c(mu.n, sigmasq.n))
# }
# 
# p.den <- function(mu.g, sigmasq.g, mu.n, sigmasq.n, Y, mu0, sigmaSq0, a, b){
#   term.1 <- mu.update.den(mu.n, Y, sigmasq.g, mu0, sigmaSq0)
#   term.2 <- sigmaSq.update.den(sigmasq.n, Y, mu.n, a, b)
#   return(term.1 * term.2)
# }
# 
# #-- --
# #MCMC function
# 
# MCMC <- function(Y, mu.init.1, sigmasq.init.1, mu.init.2, sigmasq.init.2,
#                  mu0, sigmasq0, a, b, iters){
#   
#   #chain initiation
#   mu.1 <- mu.init.1
#   sigmasq.1 <- sigmasq.init.1
#   
#   mu.2 <- mu.init.2
#   sigmasq.2 <- sigmasq.init.2
#   
#   #define storage
#   store.mu.1 <- numeric(length = iters)
#   store.sigmasq.1 <- numeric(length = iters)
#   
#   store.mu.2 <- numeric(length = iters)
#   store.sigmasq.2 <- numeric(length = iters)
#   
#   #start MCMC
#   for(i in 1:iters){
#     
#     xy.1 <- p.samp(mu.1, sigmasq.1, Y, mu0, sigmasq0, a, b)
#     mu.1.prev <- mu.1
#     sigmasq.1.prev <- sigmasq.1
#     mu.1 <- xy.1[1]
#     sigmasq.1 <- xy.1[2]
#     
#     w.upper.1 <- p.den(mu.1.prev, sigmasq.1.prev, mu.1, sigmasq.1, Y, mu0, sigmasq0, a, b)
#     w.given.xy.1 <- runif(1, 0, w.upper.1)
#     
#     if.val.1 <- p.den(mu.2, sigmasq.2, mu.1, sigmasq.1, Y, mu0, sigmasq0, a, b)
#     if (w.given.xy.1 <= if.val.1){
#       mu.2 <- mu.1
#       sigmasq.2 <- sigmasq.1
#     }else{
#       accept <- 0
#       while(accept == 0){
#         xy.2 <- p.samp(mu.2, sigmasq.2, Y, mu0, sigmasq0, a, b)
#         mu.2.cand <- xy.2[1]
#         sigmasq.2.cand <- xy.2[2]
#         
#         w.upper.2 <- p.den(mu.2, sigmasq.2, mu.2.cand, sigmasq.2.cand, Y, mu0, sigmasq0, a, b)
#         w.given.xy.2 <- runif(1, 0, w.upper.2)
#         
#         if.val.2 <- p.den(mu.1.prev, sigmasq.1.prev, mu.2.cand, sigmasq.2.cand, Y, mu0, sigmasq0, a, b)
#         if (w.given.xy.2 > if.val.2){
#           mu.2 <- mu.2.cand
#           sigmasq.2 <- sigmasq.2.cand
#           accept = 1
#         }
#       }
#     }
#     
#     store.mu.1[i] <- mu.1
#     store.sigmasq.1[i] <- sigmasq.1
#     
#     store.mu.2[i] <- mu.2
#     store.sigmasq.2[i] <- sigmasq.2
#   }
#   
#   #return chains
#   out <- list(store.mu.1 = store.mu.1,
#               store.sigmasq.1 = store.sigmasq.1,
#               store.mu.2 = store.mu.2,
#               store.sigmasq.2 = store.sigmasq.2)
#   return(out)}
# 
# #-- --
# #running the chain
# 
# 
# mu.init.1 = rnorm(1)
# sigmaSq.init.1 = rexp(1)
# mu.init.2 = rnorm(1)
# sigmaSq.init.2 = rexp(1)
# 
# MCMC.out <- MCMC(Y = Y,
#                  mu.init.1 = mu.init.1,
#                  sigmasq.init.1 = sigmaSq.init.1,
#                  mu.init.2 = mu.init.2,
#                  sigmasq.init.2 = sigmaSq.init.2,
#                  mu0 = 0,
#                  sigmasq0 = 100,
#                  a = 0.1,
#                  b = 0.1,
#                  iters = 10000)
# 
# #plotting
# par(mfrow = c(1,2))
# plot(MCMC.out$store.mu.1, xlab = "Iteration", ylab = "mu", type = "l",
#      main = paste0("First chain for mu; Mean = ", round(mean(MCMC.out$store.mu.1),3)))
# abline(a = mean(MCMC.out$store.mu.1), b = 0, lwd = 2, col = "red")
# 
# plot(MCMC.out$store.mu.2, xlab = "Iteration", ylab = "mu", type = "l",
#      main = paste0("Second chain for mu; Mean = ", round(mean(MCMC.out$store.mu.2),3)))
# abline(a = mean(MCMC.out$store.mu.2), b = 0, lwd = 2, col = "red")
# 
# plot(MCMC.out$store.sigmasq.1, xlab = "Iteration", ylab = "sigmaSq", type = "l",
#      main = paste0("First chain for sigma.sq; Mean = ", round(mean(MCMC.out$store.sigmasq.1),3)))
# abline(a = mean(MCMC.out$store.sigmasq.1), b = 0, lwd = 2, col = "red")
# 
# plot(MCMC.out$store.sigmasq.1, xlab = "Iteration", ylab = "sigmaSq", type = "l",
#      main = paste0("First chain for sigma.sq; Mean = ", round(mean(MCMC.out$store.sigmasq.1),3)))
# abline(a = mean(MCMC.out$store.sigmasq.1), b = 0, lwd = 2, col = "red")
# 
# par(mfrow = c(1,1))
# 
# #estimated value of the parameters
# mean(MCMC.out$store.mu.1)
# mean(MCMC.out$store.mu.2)
# mean(MCMC.out$store.sigmasq.1)
# mean(MCMC.out$store.sigmasq.2)
# 
# #number of couplings
# mean(MCMC.out$store.mu.1 == MCMC.out$store.mu.2)
# mean(MCMC.out$store.sigmasq.1 == MCMC.out$store.sigmasq.2)     #both should be same
# 
# #plotting joint chains
# 
# library(ggplot2)
# library(tidyverse)
# library(ggExtra)
# library(gridExtra)
# library(latex2exp)
# 
# theme_set(theme_bw(12))
# 
# df.mu <- data.frame(MCMC.out$store.mu.1, MCMC.out$store.mu.2) 
# colnames(df.mu) <- c("Chain.1", "Chain.2")
# 
# p <- ggplot(df.mu, aes(x = Chain.1, y = Chain.2)) +
#   geom_point(size = 0.8, alpha = 0.5) + 
#   geom_vline( aes(xintercept = mean(MCMC.out$store.mu.1)), col = " darkorchid3", lwd  = 0.8, lty = "dashed") +
#   geom_hline( aes(yintercept = mean(MCMC.out$store.mu.2)), col = " darkorchid3", lwd  = 0.8, lty = "dashed") +
#   labs(x = TeX(r"($\mu_1$ )"), y = TeX(r"($\mu_2$ )")) +
#   theme(axis.text = element_text(size = 13),
#         axis.title = element_text(size = 16))
# 
# ggMarginal(p, type = "histogram", fill = "#9bd5e8", col = "#033e52",
#            binwidth = 0.5, size = 3, margins = "both")
# 
# df.sigma <- data.frame(MCMC.out$store.sigmasq.1, MCMC.out$store.sigmasq.2) 
# colnames(df.sigma) <- c("Chain.1", "Chain.2")
# 
# p <- ggplot(df.sigma, aes(x = Chain.1, y = Chain.2)) +
#   geom_point(size = 0.8, alpha = 0.5) +
#   geom_vline( aes(xintercept = mean(MCMC.out$store.sigmasq.1)), col = " darkorchid3", lwd  = 0.8, lty = "dashed") +
#   geom_hline( aes(yintercept = mean(MCMC.out$store.sigmasq.2)), col = " darkorchid3", lwd  = 0.8, lty = "dashed") +
#   labs(x = TeX(r"($\sigma^2_1$ )"), y = TeX(r"($\sigma^2_2$ )")) +
#   theme(axis.text = element_text(size = 13),
#         axis.title = element_text(size = 16))
# 
# ggMarginal(p, type = "histogram", fill = "#9bd5e8", col = "#033e52",
#                   binwidth = 0.5, size = 3, margins = "both")

#===============================================================================
# Multivariate normal sample draw

# target: Multivariate Normal distribution
#         Mean: rep(2,d)
#         Sigma: CAR structure

rm(list = ls())

#-- -- -- -- -- -- -- -- --
#necessary functions

sampling.mvtn <- function(mu, sigma, samp.cur.fn){

  for (i in 1:d){
    mu.1 <- mu[i]
    mu.2 <- mu[-i]

    sigma.11 <- sigma[i,i]
    sigma.22 <- sigma[-i,-i]
    sigma.12 <- sigma[i,-i]
    sigma.21 <- sigma[-i,i]

    a <- samp.cur.fn[-i]

    i.mean <- mu.1 + sigma.12 %*% solve(sigma.22) %*% (a - mu.2)
    i.var <- sigma.11 - sigma.12 %*% solve(sigma.22) %*% sigma.21

    samp.cur.fn[i] <- rnorm(1, i.mean, sqrt(i.var))

  }

  return(samp.cur.fn)

}

density.mvtn <- function(mu, sigma, samp.cur.fn, samp.cur.fn.g){

  ret <- numeric(length = length(samp.cur.fn))

  for (i in 1:d){
    mu.1 <- mu[i]
    mu.2 <- mu[-i]

    sigma.11 <- sigma[i,i]
    sigma.22 <- sigma[-i,-i]
    sigma.12 <- sigma[i,-i]
    sigma.21 <- sigma[-i,i]
    
    if (i == 1){
      a <- samp.cur.fn.g[-i]
    }else if (i == d){
      a <- samp.cur.fn[-i]
    }else{
      a <- c(samp.cur.fn[1:(i-1)], samp.cur.fn.g[(i+1):d])
    }

    i.mean <- mu.1 + sigma.12 %*% solve(sigma.22) %*% (a - mu.2)
    i.var <- sigma.11 - sigma.12 %*% solve(sigma.22) %*% sigma.21

    ret[i] <- dnorm(samp.cur.fn[i], i.mean, sqrt(i.var))

  }

  return(prod(ret))

}

#-- -- -- -- -- -- -- -- --
#Gibbs Sampler

d <- 5

target.mean <- rgamma(d, 5, 1)                         #target parameters
foo <- matrix(rep(1:d, each = d), ncol = d)
target.sigma <- 0.5 ^ abs(foo - t(foo))

B <- 1e4
store.samp.1 <- matrix(0, ncol = d, nrow = B)
store.samp.2 <- matrix(0, ncol = d, nrow = B)

samp.cur.1 <- rnorm(d)
samp.cur.2 <- rnorm(d)

count.coupling <- 0

for (b in 1:B){

  samp.prev.1 <- samp.cur.1
  samp.cur.1 <- sampling.mvtn(target.mean, target.sigma, samp.cur.1)

  w.upper.1 <- density.mvtn(target.mean, target.sigma, samp.cur.1, samp.prev.1)
  w.given.xy.1 <- runif(1, 0, w.upper.1)

  if.val.1 <- density.mvtn(target.mean, target.sigma, samp.cur.1, samp.cur.2)
  if (w.given.xy.1 <= if.val.1){
    samp.cur.2 <- samp.cur.1
  }else{
    accept <- 0
    while(accept == 0){
      samp.cand.2 <- sampling.mvtn(target.mean, target.sigma, samp.cur.2)

      w.upper.2 <- density.mvtn(target.mean, target.sigma, samp.cand.2, samp.cur.2)
      w.given.xy.2 <- runif(1, 0, w.upper.2)

      if.val.2 <- density.mvtn(target.mean, target.sigma, samp.cand.2, samp.prev.1)
      if (w.given.xy.2 > if.val.2){
        samp.cur.2 <- samp.cand.2
        accept = 1
      }
    }
  }

  #updating coupling counts
  if ( sum(samp.cur.1 != samp.cur.2) == 0 ){
    count.coupling <- count.coupling + 1
  }

  #storing the samples
  store.samp.1[b, ] <- samp.cur.1
  store.samp.2[b, ] <- samp.cur.2
}

#-- -- -- -- -- -- -- -- --
#plotting

par(mfrow = c(1,2))

for (i in 1:d){

  plot(store.samp.1[,i], type = "l",
       main = paste0("Component: ", i, "; Chain: 1"),
       xlab = "Iterations", ylab = paste0("Component ", i),
       sub = paste0("True Val = ", round(target.mean[i], 4),
                    "; Estimated: ", round(mean(store.samp.1[,i]), 4)))
  abline(a = mean(store.samp.1[,i]), b = 0, col = "red", lwd = 2)

  plot(store.samp.2[,i], type = "l",
       main = paste0("Component: ", i, "; Chain: 2"),
       xlab = "Iterations", ylab = paste0("Component ", i),
       sub = paste0("True Val = ", round(target.mean[i], 4),
                    "; Estimated: ", round(mean(store.samp.2[,i]), 4)))
  abline(a = mean(store.samp.2[,i]), b = 0, col = "red", lwd = 2)

}

par(mfrow = c(1,1))

#-- -- -- -- -- -- -- -- --
#parameter estimation

target.mean
colMeans(store.samp.1)
colMeans(store.samp.2)

par(mfrow = c(1,2))
plot(density(as.vector(cov(store.samp.1) - target.sigma)),
     ylab = "Density",
     xlab = "Errors in covariance form first chain",
     main = "")
plot(density(as.vector(cov(store.samp.2) - target.sigma)),
     ylab = "Density",
     xlab = "Errors in covariance form second chain",
     main = "")
par(mfrow = c(1,1))

#-- -- -- -- -- -- -- -- --
#proportion of couplings

count.coupling  / B
