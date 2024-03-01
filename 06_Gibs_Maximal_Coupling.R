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

Y <- rnorm(100, mean = mu, sd = sqrt(sigmaSq))

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

#-- --
#MCMC function

Y = Y
mu.init = mean(Y)
sigmaSq.init = var(Y)
mu0 = 0
sigmaSq0 = 100
a = 0.1
b = 0.1
iters = 30000

MCMC <- function(Y, mu.init, sigmaSq.init, mu0, sigmaSq0, a, b, iters){
  
  #chain initiation
  mu <- mu.init
  sigmaSq <- sigmaSq.init
  
  #define storage
  store.mu <- numeric(length = iters)
  store.sigmasq <- numeric(length = iters)
  
  #start MCMC
  for(i in 1:iters){
    
    mu <- mu.update.samp(Y, sigmaSq, mu0, sigmaSq0)
    w.given.x <- runif(1, 0, mu.update.den(mu, Y, sigmaSq, mu0, sigmaSq0) )
    
    if (w.given.x <= sigmaSq.update.den(mu, Y, mu, a, b)){
      sigmaSq <- mu
    }else{
      accept <- 0
      while(accept == 0){
        sigmaSq.cand <- sigmaSq.update.samp(Y, mu, a, b)
        w.given.sigmaSq.cand <- runif(1, 0, sigmaSq.update.den(sigmaSq.cand, Y, mu, a, b))
        
        if (w.given.sigmaSq.cand > mu.update.den(sigmaSq.cand, Y, sigmaSq, mu0, sigmaSq0)){
          sigmaSq <- sigmaSq.cand
          accept = 1
        }
      }
    }
    
    store.mu[i] <- mu
    store.sigmasq[i] <- sigmaSq
  }
  
  #return chains
  out <- list(store.mu = store.mu,
              store.sigmasq = store.sigmasq)
  return(out)}

#-- --
#running the chain

MCMC.out <- MCMC(Y = Y,
                 mu.init = mean(Y),
                 sigmaSq.init = var(Y),
                 mu0 = 0,
                 sigmaSq0 = 100,
                 a = 0.1,
                 b = 0.1,
                 iters = 30000)

#plotting
plot(MCMC.out$store.mu, xlab = "Iteration", ylab = "mu", type = "l",
     main = paste0("Chain for mu; Mean = ", round(mean(MCMC.out$store.mu),3)))
abline(a = mean(MCMC.out$store.mu), b = 0, lwd = 2, col = "red")

plot(MCMC.out$store.sigmasq, xlab = "Iteration", ylab = "sigmaSq", type = "l",
     main = paste0("Chain for sigma.sq; Mean = ", round(mean(MCMC.out$store.sigmasq),3)))
abline(a = mean(MCMC.out$store.sigmasq), b = 0, lwd = 2, col = "red")

#estimated value of the parameters
mean(MCMC.out$store.mu)
mean(MCMC.out$store.sigmasq)

#no of couplings
mean(MCMC.out$store.mu == MCMC.out$store.sigmasq)
