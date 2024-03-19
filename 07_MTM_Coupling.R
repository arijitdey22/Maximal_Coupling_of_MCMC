#Target: Normal distribution with parameters 2 and 5

rm(list = ls())

set.seed(100)

#-- -- --
#necessary functions

q.x.y.dens <- function(x, y){
  dnorm(y, x, cand.sd)
}

q.x..samp <- function(x, n = 1){
  rnorm(n, x, cand.sd)
}

w.x.y <- function(x, y){
  lamb.x.y <- 2 / ( q.x.y.dens(x,y) + q.x.y.dens(y,x) )
  pi.x <- dnorm(x,2,sqrt(5))
  pi.x * q.x.y.dens(x,y) * lamb.x.y
}

B <- 2e4
burn.in <- 4e3
k <- 10

cand.sd <- 1

store.x.1 <- numeric(length = B)
store.x.2 <- numeric(length = B)

count.couple <- 0
count.couple.proposal <- 0

count.accept.1 <- 0
count.accept.2 <- 0

x.1 <- rnorm(1)
x.2 <- rnorm(1)

for (i in 1:B){

  coupled.proposal <- 0

  #proposing x.1.prime
  y.all.1 <- q.x..samp(x.1, k)
  w.y.all.x.1 <- w.x.y(y.all.1, x.1)

  x.1.prime <- sample(y.all.1, 1, prob = w.y.all.x.1)

  #proposing x.2.prime
  U <- runif(1)
  if (U * q.x.y.dens(x.1, x.1.prime) <= q.x.y.dens(x.2, x.1.prime)){
    x.2.prime <- x.1.prime
    coupled.proposal <- 1
  }else{

    accept = 0

    while (accept == 0){

      y.all.2 <- q.x..samp(x.2, k)
      w.y.all.x.2 <- w.x.y(y.all.2, x.2)

      x.2.tilde <- sample(y.all.2, 1, prob = w.y.all.x.2)
      V <- runif(1)

      if (V * q.x.y.dens(x.2, x.2.tilde) > q.x.y.dens(x.1, x.2.tilde)){
        accept = 1
        x.2.prime <- x.2.tilde
      }

    }

  }

  #updating the first chain
  x.star.all.1 <- q.x..samp(x.1.prime, k-1)
  x.star.all.1 <- c(x.star.all.1, x.1)

  num.1 <- sum(w.x.y(y.all.1, x.1))
  den.1 <- sum(w.x.y(x.star.all.1, x.1.prime))
  r.g.1 <- min(1, num.1 / den.1)

  B.1 <- rbinom(1, 1, r.g.1)
  x.1 <- x.1.prime * B.1 + x.1 * (1-B.1)

  #updating the second chain
  x.star.all.2 <- q.x..samp(x.2.prime, k-1)
  x.star.all.2 <- c(x.star.all.2, x.2)

  if (coupled.proposal == 1){
    y.all.2 <- q.x..samp(x.2, k)
  }
  num.2 <- sum(w.x.y(y.all.2, x.2))
  den.2 <- sum(w.x.y(x.star.all.2, x.2.prime))
  r.g.2 <- min(1, num.2 / den.2)

  B.2 <- rbinom(1, 1, r.g.2)
  x.2 <- x.2.prime * B.2 + x.2 * (1-B.2)

  #counting acceptance
  count.accept.1 <- count.accept.1 + B.1
  count.accept.2 <- count.accept.2 + B.2

  #counting coupling
  count.couple.proposal <- count.couple.proposal + coupled.proposal
  if (x.1 == x.2){
    count.couple <- count.couple + 1
  }

  #storing the chains
  store.x.1[i] <- x.1
  store.x.2[i] <- x.2

  #aesthetics
  if (i %% 1000 == 0){
    print(paste0(i," iterations completed."))
  }
}

store.x.1 <- store.x.1[ -(1:burn.in)]
store.x.2 <- store.x.2[ -(1:burn.in)]

#-- -- --
#plotting

par(mfrow = c(1,2))

plot(store.x.1, type = "l",
     main = paste0("MCMC first chain; Mean = ",
                   round(mean(store.x.1),3)),
     xlab = "Iterations",
     ylab = "X")
abline(a = mean(store.x.1), b = 0, col = "red", lwd = 2)

plot(store.x.2, type = "l",
     main = paste0("MCMC first chain; Mean = ",
                   round(mean(store.x.2),3)),
     xlab = "Iterations",
     ylab = "X")
abline(a = mean(store.x.2), b = 0, col = "red", lwd = 2)

par(mfrow = c(1,1))

#-- -- --
#parameter estimation
mean(store.x.1)
mean(store.x.2)

var(store.x.1)
var(store.x.2)

#-- -- --
#coupling proportion

count.couple.proposal / B
count.couple / B

#-- -- --
#acceptance ratio
count.accept.1 / B
count.accept.2 / B


#==============================================================================
#Target: Multivariate Normal distribution with
#        mu <- rep(2,d)
#        sigma <- 5* diag(d)

rm(list = ls())

set.seed(100)

#-- -- --
#necessary functions

q.x.y.dens <- function(x, y){
  prod(dnorm(y, x, cand.sd))
}

q.x..samp <- function(x, n = 1){                #columns represent one sample
  t(replicate(n, rnorm(d, x, cand.sd)))
}

w.x.y <- function(x, y){
  
  fun <- function(samp, y){
    lamb.x.y <- 2 / ( q.x.y.dens(samp,y) + q.x.y.dens(y,samp) )
    pi.x <- prod(dnorm(samp,2,sqrt(5)))
    pi.x * q.x.y.dens(samp,y) * lamb.x.y
  }
  
  apply(x, 1, FUN = fun, y = y)
  
}

d <- 5               #specifying the dimension

B <- 2e4
burn.in <- 4e3
k <- 10

cand.sd <- 1

store.x.1 <- matrix(0, ncol = d, nrow = B)
store.x.2 <- matrix(0, ncol = d, nrow = B)

count.couple <- 0
count.couple.proposal <- 0

count.accept.1 <- 0
count.accept.2 <- 0

x.1 <- rnorm(d)
x.2 <- rnorm(d)

for (i in 1:B){
  
  coupled.proposal <- 0
  
  #proposing x.1.prime
  y.all.1 <- q.x..samp(x.1, k)
  w.y.all.x.1 <- w.x.y(y.all.1, x.1)
  
  x.1.prime <- y.all.1[ sample(1:nrow(y.all.1), 1, prob = w.y.all.x.1), ]
  
  #proposing x.2.prime
  U <- runif(1)
  if (U * q.x.y.dens(x.1, x.1.prime) <= q.x.y.dens(x.2, x.1.prime)){
    x.2.prime <- x.1.prime
    coupled.proposal <- 1
  }else{
    
    accept = 0
    
    while (accept == 0){
      
      y.all.2 <- q.x..samp(x.2, k)
      w.y.all.x.2 <- w.x.y(y.all.2, x.2)
      
      x.2.tilde <- y.all.2[ sample(1:nrow(y.all.2), 1, prob = w.y.all.x.2), ]
      V <- runif(1)
      
      if (V * q.x.y.dens(x.2, x.2.tilde) > q.x.y.dens(x.1, x.2.tilde)){
        accept = 1
        x.2.prime <- x.2.tilde
      }
      
    }
    
  }
  
  #updating the first chain
  x.star.all.1 <- q.x..samp(x.1.prime, k-1)
  x.star.all.1 <- rbind(x.star.all.1, x.1)
  
  num.1 <- sum(w.x.y(y.all.1, x.1))
  den.1 <- sum(w.x.y(x.star.all.1, x.1.prime))
  r.g.1 <- min(1, num.1 / den.1)
  
  B.1 <- rbinom(1, 1, r.g.1)
  x.1 <- x.1.prime * B.1 + x.1 * (1-B.1)
  
  #updating the second chain
  x.star.all.2 <- q.x..samp(x.2.prime, k-1)
  x.star.all.2 <- rbind(x.star.all.2, x.2)
  
  if (coupled.proposal == 1){
    y.all.2 <- q.x..samp(x.2, k)
  }
  num.2 <- sum(w.x.y(y.all.2, x.2))
  den.2 <- sum(w.x.y(x.star.all.2, x.2.prime))
  r.g.2 <- min(1, num.2 / den.2)
  
  B.2 <- rbinom(1, 1, r.g.2)
  x.2 <- x.2.prime * B.2 + x.2 * (1-B.2)
  
  #counting acceptance
  count.accept.1 <- count.accept.1 + B.1
  count.accept.2 <- count.accept.2 + B.2
  
  #counting coupling
  count.couple.proposal <- count.couple.proposal + coupled.proposal
  if (sum(x.1 != x.2) == 0){
    count.couple <- count.couple + 1
  }
  
  #storing the chains
  store.x.1[i, ] <- x.1
  store.x.2[i, ] <- x.2
  
  #aesthetics
  if (i %% 1000 == 0){
    print(paste0(i," iterations completed."))
  }
}

store.x.1 <- store.x.1[ -(1:burn.in), ]
store.x.2 <- store.x.2[ -(1:burn.in), ]

#-- -- --
#plotting

par(mfrow = c(1,2))

for (i in 1:d){
  
  plot(store.x.1[,i], type = "l",
       main = paste0("MCMC first chain; Mean = ",
                     round(mean(store.x.1[,i]),3)),
       xlab = "Iterations",
       ylab = "X",
       sub = paste0("Component: ", i))
  abline(a = mean(store.x.1[,i]), b = 0, col = "red", lwd = 2)
  
  plot(store.x.2[,i], type = "l",
       main = paste0("MCMC first chain; Mean = ",
                     round(mean(store.x.2[,i]),3)),
       xlab = "Iterations",
       ylab = "X",
       sub = paste0("Component: ", i))
  abline(a = mean(store.x.2[,i]), b = 0, col = "red", lwd = 2)
  
}

par(mfrow = c(1,1))

#-- -- --
#parameter estimation
colMeans(store.x.1)
colMeans(store.x.2)

#-- -- --
#coupling proportion

count.couple.proposal / B
count.couple / B

#-- -- --
#acceptance ratio
count.accept.1 / B
count.accept.2 / B

