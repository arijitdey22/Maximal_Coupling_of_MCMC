library(scales)
library(mvtnorm)

###=============================================================================
###Comparison between meeting time of Algorithm 2 and Algorithm 3

#Target: Multivariate Normal with mean = rep(2,d) and sigma:
    #Case I: 2 * diag(d)
    #Case II: CAR structure. rho = 0.5

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

d.vec <- c(2,3,4,5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30)

reps <- 100                                              #number of replications for stabilizing results

store.algo.2 <- numeric(length = length(d.vec))
store.algo.3 <- numeric(length = length(d.vec))
store.algo.3.while.loop <- matrix(0, ncol = length(d.vec), nrow = reps)

for (i.d in 1:length(d.vec)){

  d <- d.vec[i.d]                                       #dimension

  target.mean <- rep(2,d)                               #target distribution

  # target.sigma <- 2 * diag(d)                         #case 1
  foo <- matrix(rep(1:d, each = d), ncol = d)         #case 2
  target.sigma <- 0.5 ^ abs(foo - t(foo))
  # foo <- matrix(0.9,  ncol = d, nrow = d)               #case 3
  # diag(foo) <- rep(1, d)
  # target.sigma <- foo

  target.sigma.inv <- solve(target.sigma)
  cons <- 1 / sqrt( det(2*pi *target.sigma) )

  ##--------------------------------------
  ##algorithm 02:

  cand.sd <- 1                                    #standard deviation of the proposal distribution

  count.iter.2.vec <- numeric(length = reps)      #storage for counts

  for (i.reps in 1:reps){                         #loop for reps

    #set.seed(i.reps)

    x.cur.2 <- runif(d)                           #initial values
    y.cur.2 <- runif(d)

    y.accept.alg.2 <- 0                           #marker variable for acceptance
    count.iter.2 <- 0                             #count of the first meeting

    while (y.accept.alg.2 == 0){

      count.iter.2 <- count.iter.2 + 1

      #proposal of x
      x.cand.2 <- x.cur.2 + cand.sd * rnorm(d)
      U <- runif(1)
      log.val.U <- sum(dnorm(x.cand.2, y.cur.2, cand.sd, log = T)) - sum(dnorm(x.cand.2, x.cur.2, cand.sd, log = T))

      #proposal of y
      if (log(U) <= log.val.U){
        y.cand.2 <- x.cand.2
      }else{
        accept.2 = 0
        while(accept.2 == 0){
          y.curl.2 <- y.cur.2 + cand.sd * rnorm(d)
          V <- runif(1)
          log.V.val <- sum(dnorm(y.curl.2, x.cur.2, cand.sd, log = T)) - sum(dnorm(y.curl.2, y.cur.2, cand.sd, log = T))

          if(log(V) > log.V.val){
            y.cand.2 <- y.curl.2
            accept.2 = 1
          }
        }
      }

      #acceptance probabilities
      a.xx <- fun.a.x.y(x.cur.2, x.cand.2, target.mean, target.sigma.inv, cons)
      a.yy <- fun.a.x.y(y.cur.2, y.cand.2, target.mean, target.sigma.inv, cons)

      B.x <- rbinom(1, 1, min(1, a.xx))
      B.y <- rbinom(1, 1, min(1, a.yy))

      if ( ( sum(x.cand.2 != y.cand.2) == 0 ) & (B.x == B.y) ){
        y.accept.alg.2 <- 1
      }

    }

    #aesthetics
    print(paste0("Algorithm 2; d = ", d.vec[i.d], ", rep = ", i.reps))
    print("==================================================")

    count.iter.2.vec[i.reps] <- count.iter.2
  }

  ##--------------------------------------
  ##algorithm 03:

  cand.sd <- 1                                    #standard deviation of the proposal distribution

  count.iter.3.vec <- numeric(length = reps)      #storage for counts

  for (i.reps in 1:reps){                         #loop for reps

    while.loop.count <- 0

    #set.seed(i.reps)

    x.cur.3 <- runif(d)                           #initial values
    y.cur.3 <- runif(d)

    y.accept.alg.3 <- 0                           #marker variable for acceptance
    count.iter.3 <- 0                             #count of the first meeting

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

    #aesthetics
    print(paste0("While loop count: ", while.loop.count))
    print("--------------------")
    print(paste0("Algorithm 3; d = ", d.vec[i.d], ", rep = ", i.reps))
    print("==================================================")

    count.iter.3.vec[i.reps] <- count.iter.3
    store.algo.3.while.loop[i.reps, i.d] <- while.loop.count
  }

  store.algo.2[i.d] <- mean(count.iter.2.vec)
  store.algo.3[i.d] <- mean(count.iter.3.vec)

}

##--------------------------------------

# save(store.algo.2, store.algo.3, store.algo.3.while.loop, file = "05_Comparison.i.1.Rdata")
save(store.algo.2, store.algo.3, store.algo.3.while.loop, file = "05_Comparison.c.1.Rdata")

# load("05_Comparison.i.1.Rdata")
load("05_Comparison.c.1.Rdata")

d.vec <- c(2,3,4,5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30)

library(ggplot2)

#plotting meeting times

ggplot() +
  geom_line(aes(x = d.vec, y = store.algo.2), col = "darkslategray4", lwd = 1)  +
  labs(x = "Dimensions",
       y = "Meeting time",
       col = "") + # ylim(0,max(c(store.algo.2, store.algo.3))) +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 14))

ggplot() +
  geom_line(aes(x = d.vec, y = store.algo.3), col = "darkslategray4", lwd = 1)  +
  labs(x = "Dimensions",
       y = "Meeting time",
       col = "") + ylim(0,max(c(store.algo.2, store.algo.3))) +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 14))

#plotting while loop count

algo.3.while.loop <- colMeans(store.algo.3.while.loop)

ggplot() +
  geom_line(aes(x = d.vec, y = algo.3.while.loop), col = "darkslategray4", lwd = 1)  +
  labs(x = "Dimensions",
       y = "While loop count",
       col = "") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 14))


###=============================================================================
###plot of d vs Gibbs maximal coupling time

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

d.vec <- c(2,3,4,5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30)
store.count.d <- numeric(length = length(d.vec))

for ( i.d in 1:length(d.vec)){

  d <- d.vec[i.d]

  target.mean <- rep(2,d)                               #target distribution

  target.sigma <- 2 * diag(d)                         #case 1
  # foo <- matrix(rep(1:d, each = d), ncol = d)         #case 2
  # target.sigma <- 0.5 ^ abs(foo - t(foo))
  # foo <- matrix(0.9,  ncol = d, nrow = d)             #case 3
  # diag(foo) <- rep(1, d)
  # target.sigma <- foo

  reps <- 100
  count.iter.vec <- numeric(length = reps)

  for (i.reps in 1:reps){

    samp.cur.1 <- rnorm(d)
    samp.cur.2 <- rnorm(d)

    coupling.index <- 0
    count.iter <- 0

    while (coupling.index == 0){

      count.iter <- count.iter + 1

      samp.prev.1 <- samp.cur.1
      samp.cur.1 <- sampling.mvtn(target.mean, target.sigma, samp.cur.1)

      w.upper.1 <- density.mvtn(target.mean, target.sigma, samp.cur.1, samp.prev.1)
      w.given.xy.1 <- runif(1, 0, w.upper.1)

      if.val.1 <- density.mvtn(target.mean, target.sigma, samp.cur.1, samp.cur.2)
      if (w.given.xy.1 <= if.val.1){
        samp.cur.2 <- samp.cur.1
        coupling.index <- 1
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
    }

    #aesthetics
    print(paste0("We are at d = ", d.vec[i.d],", and rep = ", i.reps))

    count.iter.vec[i.reps] <- count.iter

  }

  store.count.d[i.d] <- mean(count.iter.vec)

}

save(d.vec, store.count.d, file = "05_Comparison.i.2.Rdata")
# save(d.vec, store.count.d, file = "05_Comparison.c.2.Rdata")
#
# ##--------------------------------------

load("05_Comparison.i.2.Rdata")
# load("05_Comparison.c.2.Rdata")

library(ggplot2)

ggplot( ) +
  geom_line(aes(x = d.vec, y = store.count.d), col = "darkslategray4", lwd = 1) +
  labs(x = "Dimensions",
       y = "Meeting time",
       col = "") +
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 14))
