#Target: Normal distribution with parameters 2 and 5

#-- -- --
#necessary functions

rm(list = ls())

T.x.y.dens <- function(x, y, cand.sd){
  dnorm(y, x, cand.sd)
}

T.x..samp <- function(x, cand.sd, n = 1){
  rnorm(n, x, cand.sd)
}

w.x.y <- function(x, y, cand.sd){
  lamb.x.y <- 2 / ( T.x.y.dens(x,y,cand.sd) + T.x.y.dens(y,x,cand.sd) )
  dnorm(x,2,sqrt(5)) * T.x.y.dens(x,y,cand.sd) * lamb.x.y
}

w.star.x.y <- function(x, y, y.star, cand.sd, lam){
  lam * w.x.y(x,y,cand.sd) + (1-lam) * T.x.y.dens(x,y.star,cand.sd)
}

B <- 5e4
k <- 10
lam <- 0.5

cand.sd <- 1

store.prob.1 <- numeric(length = B)
store.prob.2 <- numeric(length = B)

x.cur.1 <- rnorm(1)
x.cur.2 <- rnorm(1)


for (i in 1:B){

  set.seed(i)
  
  #updating x

  y.all <- T.x..samp(x.cur.1, cand.sd, k)
  
  w.y.all.x <- w.x.y(y.all, x.cur.1, cand.sd)
  w.y.all.x.star <- w.star.x.y(y.all, x.cur.1, x.cur.2, cand.sd, lam)

  y <- sample(y.all, 1, prob = w.y.all.x)
  y.star <- sample(y.all, 1, prob = w.y.all.x.star)

  store.prob.1[i] <- T.x.y.dens(x.cur.2, y, cand.sd) / T.x.y.dens(x.cur.1, y, cand.sd)
  store.prob.2[i] <- T.x.y.dens(x.cur.2, y.star, cand.sd) / T.x.y.dens(x.cur.1, y.star, cand.sd)
  
  if (i %% 1000 == 0){
    print(paste0(i," iterations completed."))
  }
}


#-- -- --
#plotting

par(mfrow = c(1,2))

plot(density(store.prob.1),
     main = "Density with old algo",
     xlab = "Probability",
     ylab = "Density")
plot(density(store.prob.2),
     main = "Density with new algo",
     xlab = "Probability",
     ylab = "Density")
par(mfrow = c(1,1))

#-- -- --
#average prob
mean(store.prob.1)
mean(store.prob.2)
