rm(list = ls())

#-- -- 
#target: X ~ N(2,1)
#        Y ~ N(-1, 4)

B <- 1e5

store.x <- numeric(length = B)
store.y <- numeric(length = B)

for (i in 1:B){
  
  store.x[i] <- rnorm(1, 2, 1)
  
  w.given.x <- runif(1, 0, dnorm(store.x[i], 2, 1))
  
  if (w.given.x <= dnorm(store.x[i], -1, 2)){
    store.y[i] <- store.x[i]
  }else{
    accept <- 0
    while(accept == 0){
      y.cand <- rnorm(1, -1, 2)
      w.given.y.cand <- runif(1, 0, dnorm(y.cand, -1, 2))
      
      if (w.given.y.cand > dnorm(y.cand, 2,1)){
        store.y[i] <- y.cand
        accept = 1
      }
    }
  }
  
}

#-- -- 
#plotting

par(mfrow = c(1,2))

plot(store.x, type = "l", 
     main = paste0("MCMC chain for X, Mean = ", round(mean(store.x),3)),
     xlab = "Iterations",
     ylab = "X")
abline(a = mean(store.x), b = 0, col = "red", lwd = 2)

plot(store.y, type = "l", 
     main = paste0("MCMC chain for X, Mean = ", round(mean(store.y),3)),
     xlab = "Iterations",
     ylab = "X")
abline(a = mean(store.y), b = 0, col = "red", lwd = 2)

par(mfrow = c(1,1))

#-- -- 
#parameter estimation
mean(store.x)
var(store.x)

mean(store.y)
var(store.y)

#-- -- 
#coupling ratio

mean(store.x == store.y)
