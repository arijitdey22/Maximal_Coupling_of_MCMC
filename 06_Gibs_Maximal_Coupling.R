# rm(list = ls())
# 
# set.seed(100)
# 
# #-- -- 
# #data generations
# 
# mu <- 5
# sigmaSq <- 1
# 
# Y <- rnorm(100, mean = mu, sd = sqrt(sigmaSq))
# 
# #-- -- 
# #full conditional distributions
# #of mu given sigma
# 
# mu.update <- function(Y, sigmaSq, mu0, sigmaSq0){
#   mu.posvar <- 1 / (length(Y) / sigmaSq + 1 / sigmaSq0)
#   mu.posmean <- (sum(Y) / sigmaSq + mu0 / sigmaSq0) * mu.posvar
#   out <- mu.posmean + sqrt(mu.posvar) * rnorm(1)
#   return(out)
#   }
# 
# #of sigma given mu
# 
# sigmaSq.update <- function(Y, mu, a, b){
#   sigmaSq.inv <- rgamma(1, shape = a + length(Y) / 2,
#                         rate = b + sum((Y- mu)^2) / 2)
#   out <- 1 / sigmaSq.inv
#   return(out)
#   }
# 
# #-- --
# #MCMC function
# 
# MCMC <- function(Y, mu.init, sigmaSq.init, mu0, sigmaSq0, a, b, iters){
#   
#   #chain initiation
#   mu <- mu.init
#   sigmaSq <- sigmaSq.init
#   
#   #define storage
#   store.mu <- numeric(length = iters)
#   store.sigmasq <- numeric(length = iters)
#   
#   #start MCMC
#   for(i in 1:iters){
#     mu <- mu.update(Y, sigmaSq, mu0, sigmaSq0)
#     sigmaSq <- sigmaSq.update(Y, mu, a, b)
#     store.mu[i] <- mu
#     store.sigmasq[i] <- sigmaSq
#   }
#   
#   #return chains
#   out <- list(store.mu = store.mu,
#               store.sigmasq = store.sigmasq)
#   return(out)}
# 
# #-- --
# #running the chain
# 
# MCMC.out <- MCMC(Y = Y,
#                  mu.init = mean(Y),
#                  sigmaSq.init = var(Y),
#                  mu0 = 0,
#                  sigmaSq0 = 100,
#                  a = 0.1,
#                  b = 0.1,
#                  iters = 30000)
# 
# #plotting
# plot(MCMC.out$store.mu, xlab = "Iteration", ylab = "mu", type = "l",
#      main = paste0("Chain for mu; Mean = ", round(mean(MCMC.out$store.mu),3)))
# abline(a = mean(MCMC.out$store.mu), b = 0, lwd = 2, col = "red")
# 
# plot(MCMC.out$store.sigmasq, xlab = "Iteration", ylab = "sigmaSq", type = "l",
#      main = paste0("Chain for sigma.sq; Mean = ", round(mean(MCMC.out$store.sigmasq),3)))
# abline(a = mean(MCMC.out$store.sigmasq), b = 0, lwd = 2, col = "red")
# 
# #estimated value of the parameters
# mean(MCMC.out$store.mu)
# mean(MCMC.out$store.sigmasq)


###============================================================================
###Using maximal coupling

rm(list = ls())

library(invgamma)
set.seed(100)

#-- -- 
#data generations

mu <- 5
sigmaSq <- 1

Y <- rnorm(500, mean = mu, sd = sqrt(sigmaSq))

#-- -- 
#full conditional distributions
#of mu given sigma

mu.update.samp <- function(Y, sigmaSq, mu0, sigmaSq0){
  mu.posvar <- 1 / (length(Y) / sigmaSq + 1 / sigmaSq0)
  mu.posmean <- (sum(Y) / sigmaSq + mu0 / sigmaSq0) * mu.posvar
  out <- mu.posmean + sqrt(mu.posvar) * rnorm(1)
  return(out)
}

mu.update.den <- function(mu.fn, Y, sigmaSq, mu0, sigmaSq0){
  mu.posvar <- 1 / (length(Y) / sigmaSq + 1 / sigmaSq0)
  mu.posmean <- (sum(Y) / sigmaSq + mu0 / sigmaSq0) * mu.posvar
  out <- dnorm(mu.fn, mu.posmean, sqrt(mu.posvar))
  return(out)
}

#of sigma given mu

sigmaSq.update.samp <- function(Y, mu, a, b){
  sigmaSq.inv <- rgamma(1, shape = a + length(Y) / 2,
                        rate = b + sum((Y- mu)^2) / 2)
  out <- 1 / sigmaSq.inv
  return(out)
}

sigmaSq.update.den <- function(sigmaSq.fn, Y, mu, a, b){
  if (sigmaSq.fn <= 0){
    out <- 0
  }else{
    out <- dinvgamma(sigmaSq.fn, shape = a + length(Y) / 2,
                          rate = b + sum((Y- mu)^2) / 2)
  }

  return(out)
}

#of the p kernel

p.samp <- function(mu.g, sigmasq.g, Y, mu0, sigmaSq0, a, b){
  mu.n <- mu.update.samp(Y, sigmasq.g, mu0, sigmaSq0)
  sigmasq.n <- sigmaSq.update.samp(Y, mu.n, a, b)
  return(c(mu.n, sigmasq.n))
}

p.den <- function(mu.g, sigmasq.g, mu.n, sigmasq.n, Y, mu0, sigmaSq0, a, b){
  term.1 <- mu.update.den(mu.n, Y, sigmasq.g, mu0, sigmaSq0)
  term.2 <- sigmaSq.update.den(sigmasq.n, Y, mu.n, a, b)
  return(term.1 * term.2)
}

#-- --
#MCMC function

MCMC <- function(Y, mu.init.1, sigmasq.init.1, mu.init.2, sigmasq.init.2,
                 mu0, sigmasq0, a, b, iters){
  
  #chain initiation
  mu.1 <- mu.init.1
  sigmasq.1 <- sigmasq.init.1
  
  mu.2 <- mu.init.2
  sigmasq.2 <- sigmasq.init.2
  
  #define storage
  store.mu.1 <- numeric(length = iters)
  store.sigmasq.1 <- numeric(length = iters)
  
  store.mu.2 <- numeric(length = iters)
  store.sigmasq.2 <- numeric(length = iters)
  
  #start MCMC
  for(i in 1:iters){
    
    xy.1 <- p.samp(mu.1, sigmasq.1, Y, mu0, sigmasq0, a, b)
    mu.1.prev <- mu.1
    sigmasq.1.prev <- sigmasq.1
    mu.1 <- xy.1[1]
    sigmasq.1 <- xy.1[2]
    
    w.upper.1 <- p.den(mu.1.prev, sigmasq.1.prev, mu.1, sigmasq.1, Y, mu0, sigmasq0, a, b)
    w.given.xy.1 <- runif(1, 0, w.upper.1)
    
    if.val.1 <- p.den(mu.2, sigmasq.2, mu.1, sigmasq.1, Y, mu0, sigmasq0, a, b)
    if (w.given.xy.1 <= if.val.1){
      mu.2 <- mu.1
      sigmasq.2 <- sigmasq.1
    }else{
      accept <- 0
      while(accept == 0){
        xy.2 <- p.samp(mu.2, sigmasq.2, Y, mu0, sigmasq0, a, b)
        mu.2.cand <- xy.2[1]
        sigmasq.2.cand <- xy.2[2]
        
        w.upper.2 <- p.den(mu.2, sigmasq.2, mu.2.cand, sigmasq.2.cand, Y, mu0, sigmasq0, a, b)
        w.given.xy.2 <- runif(1, 0, w.upper.2)
        
        if.val.2 <- p.den(mu.1.prev, sigmasq.1.prev, mu.2.cand, sigmasq.2.cand, Y, mu0, sigmasq0, a, b)
        if (w.given.xy.2 > if.val.2){
          mu.2 <- mu.2.cand
          sigmasq.2 <- sigmasq.2.cand
          accept = 1
        }
      }
    }
    
    store.mu.1[i] <- mu.1
    store.sigmasq.1[i] <- sigmasq.1
    
    store.mu.2[i] <- mu.2
    store.sigmasq.2[i] <- sigmasq.2
  }
  
  #return chains
  out <- list(store.mu.1 = store.mu.1,
              store.sigmasq.1 = store.sigmasq.1,
              store.mu.2 = store.mu.2,
              store.sigmasq.2 = store.sigmasq.2)
  return(out)}

#-- --
#running the chain


mu.init.1 = rnorm(1)
sigmaSq.init.1 = rexp(1)
mu.init.2 = rnorm(1)
sigmaSq.init.2 = rexp(1)

MCMC.out <- MCMC(Y = Y,
                 mu.init.1 = mu.init.1,
                 sigmasq.init.1 = sigmaSq.init.1,
                 mu.init.2 = mu.init.2,
                 sigmasq.init.2 = sigmaSq.init.2,
                 mu0 = 0,
                 sigmasq0 = 100,
                 a = 0.1,
                 b = 0.1,
                 iters = 30000)

#plotting
par(mfrow = c(1,2))
plot(MCMC.out$store.mu.1, xlab = "Iteration", ylab = "mu", type = "l",
     main = paste0("First chain for mu; Mean = ", round(mean(MCMC.out$store.mu.1),3)))
abline(a = mean(MCMC.out$store.mu.1), b = 0, lwd = 2, col = "red")

plot(MCMC.out$store.mu.2, xlab = "Iteration", ylab = "mu", type = "l",
     main = paste0("Second chain for mu; Mean = ", round(mean(MCMC.out$store.mu.2),3)))
abline(a = mean(MCMC.out$store.mu.2), b = 0, lwd = 2, col = "red")

plot(MCMC.out$store.sigmasq.1, xlab = "Iteration", ylab = "sigmaSq", type = "l",
     main = paste0("First chain for sigma.sq; Mean = ", round(mean(MCMC.out$store.sigmasq.1),3)))
abline(a = mean(MCMC.out$store.sigmasq.1), b = 0, lwd = 2, col = "red")

plot(MCMC.out$store.sigmasq.1, xlab = "Iteration", ylab = "sigmaSq", type = "l",
     main = paste0("First chain for sigma.sq; Mean = ", round(mean(MCMC.out$store.sigmasq.1),3)))
abline(a = mean(MCMC.out$store.sigmasq.1), b = 0, lwd = 2, col = "red")

par(mfrow = c(1,1))

#estimated value of the parameters
mean(MCMC.out$store.mu.1)
mean(MCMC.out$store.mu.2)
mean(MCMC.out$store.sigmasq.1)
mean(MCMC.out$store.sigmasq.2)

#number of couplings
mean(MCMC.out$store.mu.1 == MCMC.out$store.mu.2)
mean(MCMC.out$store.sigmasq.1 == MCMC.out$store.sigmasq.2)     #both should be same

