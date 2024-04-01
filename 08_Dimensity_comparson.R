library(scales)
library(mvtnorm)

###=============================================================================
###Comparison between meeting time of Algorithm 2, Algorithm 3, 
###                and MTM coupling Algorithm.

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

q.x.y.dens <- function(x, y, cand.sd){
  prod(dnorm(y, x, cand.sd))
}

q.x..samp <- function(x, n = 1, d, cand.sd){                #columns represent one sample
  t(replicate(n, rnorm(d, x, cand.sd)))
}

w.x.y <- function(x, y, mu, sigma.inv, cons, cand.sd){
  
  fun <- function(samp, y, mu, sigma.inv, cons, cand.sd){
    lamb.x.y <- 2 / ( q.x.y.dens(samp, y, cand.sd) + q.x.y.dens(y, samp, cand.sd) )
    pi.x <- my.dmvnorm(samp, mu, sigma.inv, cons)
    pi.x * q.x.y.dens(samp, y, cand.sd) * lamb.x.y
  }
  
  apply(x, 1, FUN = fun, y = y, mu = mu, sigma.inv = sigma.inv, cons = cons, cand.sd = cand.sd)
  
}

##--------------------------------------

d.vec <- c(2,3,4,5, 6, 7, 8, 9, 10, 12, 15)

reps <- 100                                              #number of replications for stabilizing results

store.algo.2 <- numeric(length = length(d.vec))
store.algo.3 <- numeric(length = length(d.vec))
store.algo.MTM <- numeric(length = length(d.vec))

store.algo.2.time <- numeric(length = length(d.vec))
store.algo.3.time <- numeric(length = length(d.vec))
store.algo.MTM.time <- numeric(length = length(d.vec))

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
  time.iter.2.vec <- numeric(length = reps)       #storage for time
  
  for (i.reps in 1:reps){                         #loop for reps
    
    tick <- proc.time()[3]
    
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
    
    tock <- proc.time()[3]
    
    count.iter.2.vec[i.reps] <- count.iter.2
    time.iter.2.vec[i.reps] <- tock - tick
  }
  
  ##--------------------------------------
  ##algorithm 03:
  
  cand.sd <- 1                                    #standard deviation of the proposal distribution
  
  count.iter.3.vec <- numeric(length = reps)      #storage for counts
  time.iter.3.vec <- numeric(length = reps)       #storage for time
  
  for (i.reps in 1:reps){                         #loop for reps
    
    tick <- proc.time()[3]
    
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
    print(paste0("Algorithm 3; d = ", d.vec[i.d], ", rep = ", i.reps))
    print("==================================================")
    
    tock <- proc.time()[3]
    
    count.iter.3.vec[i.reps] <- count.iter.3
    time.iter.3.vec[i.reps] <- tock - tick
    
  }
  
  ##--------------------------------------
  ##algorithm MTM:
  
  cand.sd <- 1                                      #standard deviation of the proposal distribution
  k <- 10
  
  count.iter.MTM.vec <- numeric(length = reps)      #storage for counts
  time.iter.MTM.vec <- numeric(length = reps)       #storage for time
  
  for (i.reps in 1:reps){
    
    tick <- proc.time()[3]
    
    #set.seed(i.reps)
    
    x.1 <- runif(d)                             #initial values
    x.2 <- runif(d)
    
    y.accept.alg.MTM <- 0                           #marker variable for acceptance
    count.iter.MTM <- 0                             #count of the first meeting
    
    while (y.accept.alg.MTM == 0){
      
      count.iter.MTM <- count.iter.MTM + 1
      
      #proposing x.1.prime
      y.all.1 <- q.x..samp(x.1, k, d, cand.sd)
      w.y.all.x.1 <- w.x.y(y.all.1, x.1, target.mean, target.sigma.inv, cons, cand.sd)
      
      x.1.prime <- y.all.1[ sample(1:nrow(y.all.1), 1, prob = w.y.all.x.1), ]
      
      #proposing x.2.prime
      U <- runif(1)
      if (U * q.x.y.dens(x.1, x.1.prime, cand.sd) <= q.x.y.dens(x.2, x.1.prime, cand.sd)){
        x.2.prime <- x.1.prime
        y.all.2 <- y.all.1
      }else{
        
        accept = 0
        
        while (accept == 0){
          
          y.all.2 <- q.x..samp(x.2, k, d, cand.sd)
          w.y.all.x.2 <- w.x.y(y.all.2, x.2, target.mean, target.sigma.inv, cons, cand.sd)
          
          x.2.tilde <- y.all.2[ sample(1:nrow(y.all.2), 1, prob = w.y.all.x.2), ]
          V <- runif(1)
          
          if (V * q.x.y.dens(x.2, x.2.tilde, cand.sd) > q.x.y.dens(x.1, x.2.tilde, cand.sd)){
            accept = 1
            x.2.prime <- x.2.tilde
          }
          
        }
        
      }
      
      #common acceptance value
      U <- runif(1)
      
      #updating the first chain
      x.star.all <- q.x..samp(x.1.prime, k-1, d, cand.sd)
      x.star.all.1 <- rbind(x.star.all, x.1)
      
      num.1 <- sum(w.x.y(y.all.1, x.1, target.mean, target.sigma.inv, cons, cand.sd))
      den.1 <- sum(w.x.y(x.star.all.1, x.1.prime, target.mean, target.sigma.inv, cons, cand.sd))
      r.g.1 <- min(1, num.1 / den.1)
      
      B.1 <- r.g.1 < U
      x.1 <- x.1.prime * B.1 + x.1 * (1-B.1)
      
      #updating the second chain
      x.star.all.2 <- rbind(x.star.all, x.2)
      
      num.2 <- sum(w.x.y(y.all.2, x.2, target.mean, target.sigma.inv, cons, cand.sd))
      den.2 <- sum(w.x.y(x.star.all.2, x.2.prime, target.mean, target.sigma.inv, cons, cand.sd))
      r.g.2 <- min(1, num.2 / den.2)
      
      B.2 <- r.g.2 < U
      x.2 <- x.2.prime * B.2 + x.2 * (1-B.2)
      
      if (sum(x.1 != x.2) == 0){
        y.accept.alg.MTM <- 1
      }
      
    }
    
    #aesthetics
    print(paste0("Algorithm MTM; d = ", d.vec[i.d], ", rep = ", i.reps))
    print("==================================================")
    
    tock <- proc.time()[3]
    
    count.iter.MTM.vec[i.reps] <- count.iter.MTM
    time.iter.MTM.vec[i.reps] <- tock - tick
  
  }
  
  store.algo.2[i.d] <- mean(count.iter.2.vec)
  store.algo.3[i.d] <- mean(count.iter.3.vec)
  store.algo.MTM[i.d] <- mean(count.iter.MTM.vec)
  
  store.algo.2.time[i.d] <- mean(time.iter.2.vec)
  store.algo.3.time[i.d] <- mean(time.iter.3.vec)
  store.algo.MTM.time[i.d] <- mean(time.iter.MTM.vec)
  
}

store.algo.2
store.algo.3
store.algo.MTM
store.algo.2.time
store.algo.3.time
store.algo.MTM.time

# save(store.algo.2, store.algo.3, store.algo.MTM,
#      store.algo.2.time, store.algo.3.time, store.algo.MTM.time,
#      file = "08_Dimensity_comparson.i.Rdata")
save(store.algo.2, store.algo.3, store.algo.MTM,
     store.algo.2.time, store.algo.3.time, store.algo.MTM.time,
     file = "08_Dimensity_comparson.c.Rdata")

##--------------------------------------
## Plotting

# load("08_Dimensity_comparson.i.Rdata")
load("08_Dimensity_comparson.c.Rdata")

d.vec <- c(2,3,4,5, 6, 7, 8, 9, 10, 12, 15)

#iterations

df.1 <- data.frame(X = rep(d.vec, 3),
                   Y = c(store.algo.2, store.algo.3, store.algo.MTM),
                   type = c(rep("Algotithm 2", length(d.vec)),
                            rep("Algotithm 3", length(d.vec)),
                            rep("Algotithm MTM", length(d.vec))))

p.1 <- ggplot( data = df.1 ) +
  geom_line(aes(x = X, y = Y, col = type), lwd = 1) +
  labs(x = "Dimensions",
       y = "Meeting time (Iterations)",
       col = "") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.position = "bottom")

#-- --
#time

df.2 <- data.frame(X = rep(d.vec, 3),
                   Y = c(store.algo.2.time, store.algo.3.time, store.algo.MTM.time),
                   type = c(rep("Algotithm 2", length(d.vec)),
                            rep("Algotithm 3", length(d.vec)),
                            rep("Algotithm MTM", length(d.vec))))

p.2<- ggplot( data = df.2 ) +
  geom_line(aes(x = X, y = Y, col = type), lwd = 1) +
  labs(x = "Dimensions",
       y = "Meeting time (Real time)",
       col = "") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.position = "bottom")

library(patchwork)

p.1 + p.2 + plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

