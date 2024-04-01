#Target: Normal distribution with parameters 2 and 5

# rm(list = ls())
# 
# set.seed(100)
# 
# #-- -- --
# #necessary functions
# 
# q.x.y.dens <- function(x, y){
#   dnorm(y, x, cand.sd)
# }
# 
# q.x..samp <- function(x, n = 1){
#   rnorm(n, x, cand.sd)
# }
# 
# w.x.y <- function(x, y){
#   lamb.x.y <- 2 / ( q.x.y.dens(x,y) + q.x.y.dens(y,x) )
#   pi.x <- dnorm(x,2,sqrt(5))
#   pi.x * q.x.y.dens(x,y) * lamb.x.y
# }
# 
# B <- 2e4
# burn.in <- 0
# k <- 10
# 
# cand.sd <- 2
# 
# store.x.1 <- numeric(length = B)
# store.x.2 <- numeric(length = B)
# 
# count.couple <- 0
# count.couple.proposal <- 0
# 
# count.accept.1 <- 0
# count.accept.2 <- 0
# 
# x.1 <- rnorm(1)
# x.2 <- rnorm(1)
# 
# for (i in 1:B){
# 
#   coupled.proposal <- 0
# 
#   #proposing x.1.prime
#   y.all.1 <- q.x..samp(x.1, k)
#   w.y.all.x.1 <- w.x.y(y.all.1, x.1)
# 
#   x.1.prime <- sample(y.all.1, 1, prob = w.y.all.x.1)
# 
#   #proposing x.2.prime
#   U <- runif(1)
#   if (U * q.x.y.dens(x.1, x.1.prime) <= q.x.y.dens(x.2, x.1.prime)){
#     x.2.prime <- x.1.prime
#     y.all.2 <- y.all.1
#     coupled.proposal <- 1
#   }else{
# 
#     accept = 0
# 
#     while (accept == 0){
# 
#       y.all.2 <- q.x..samp(x.2, k)
#       w.y.all.x.2 <- w.x.y(y.all.2, x.2)
# 
#       x.2.tilde <- sample(y.all.2, 1, prob = w.y.all.x.2)
#       V <- runif(1)
# 
#       if (V * q.x.y.dens(x.2, x.2.tilde) > q.x.y.dens(x.1, x.2.tilde)){
#         accept = 1
#         x.2.prime <- x.2.tilde
#       }
# 
#     }
# 
#   }
# 
#   #common acceptance value
#   U <- runif(1)
# 
#   #updating the first chain
#   x.star.all <- q.x..samp(x.1.prime, k-1)
#   x.star.all.1 <- c(x.star.all, x.1)
# 
#   num.1 <- sum(w.x.y(y.all.1, x.1))
#   den.1 <- sum(w.x.y(x.star.all.1, x.1.prime))
#   r.g.1 <- min(1, num.1 / den.1)
# 
#   B.1 <- r.g.1 < U
#   x.1 <- x.1.prime * B.1 + x.1 * (1-B.1)
# 
#   #updating the second chain
#   x.star.all.2 <- c(x.star.all, x.2)
# 
#   num.2 <- sum(w.x.y(y.all.2, x.2))
#   den.2 <- sum(w.x.y(x.star.all.2, x.2.prime))
#   r.g.2 <- min(1, num.2 / den.2)
# 
#   B.2 <- r.g.1 < U
#   x.2 <- x.2.prime * B.2 + x.2 * (1-B.2)
# 
#   #counting acceptance
#   count.accept.1 <- count.accept.1 + as.numeric(B.1)
#   count.accept.2 <- count.accept.2 + as.numeric(B.2)
# 
#   #counting coupling
#   count.couple.proposal <- count.couple.proposal + coupled.proposal
#   if (x.1 == x.2){
#     count.couple <- count.couple + 1
#   }
# 
#   #storing the chains
#   store.x.1[i] <- x.1
#   store.x.2[i] <- x.2
# 
#   #aesthetics
#   if (i %% 1000 == 0){
#     print(paste0(i," iterations completed."))
#   }
# }
# 
# # store.x.1 <- store.x.1[ -(1:burn.in)]
# # store.x.2 <- store.x.2[ -(1:burn.in)]
# 
# #-- -- --
# #plotting
# 
# par(mfrow = c(1,2))
# 
# plot(store.x.1, type = "l",
#      main = paste0("MCMC first chain; Mean = ",
#                    round(mean(store.x.1),3)),
#      xlab = "Iterations",
#      ylab = "X")
# abline(a = mean(store.x.1), b = 0, col = "red", lwd = 2)
# 
# plot(store.x.2, type = "l",
#      main = paste0("MCMC first chain; Mean = ",
#                    round(mean(store.x.2),3)),
#      xlab = "Iterations",
#      ylab = "X")
# abline(a = mean(store.x.2), b = 0, col = "red", lwd = 2)
# 
# par(mfrow = c(1,1))
# 
# #-- -- --
# #parameter estimation
# mean(store.x.1)
# mean(store.x.2)
# 
# var(store.x.1)
# var(store.x.2)
# 
# #-- -- --
# #coupling proportion
# 
# count.couple.proposal / B
# count.couple / B
# 
# #-- -- --
# #acceptance ratio
# count.accept.1 / B
# count.accept.2 / B


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

d <- 6               #specifying the dimension

B <- 1e4
burn.in <- 0
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
    y.all.2 <- y.all.1
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

  #common acceptance value
  U <- runif(1)

  #updating the first chain
  x.star.all <- q.x..samp(x.1.prime, k-1)
  x.star.all.1 <- rbind(x.star.all, x.1)

  num.1 <- sum(w.x.y(y.all.1, x.1))
  den.1 <- sum(w.x.y(x.star.all.1, x.1.prime))
  r.g.1 <- min(1, num.1 / den.1)

  B.1 <- r.g.1 < U
  x.1 <- x.1.prime * B.1 + x.1 * (1-B.1)

  #updating the second chain
  x.star.all.2 <- rbind(x.star.all, x.2)

  num.2 <- sum(w.x.y(y.all.2, x.2))
  den.2 <- sum(w.x.y(x.star.all.2, x.2.prime))
  r.g.2 <- min(1, num.2 / den.2)

  B.2 <- r.g.2 < U
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

# store.x.1 <- store.x.1[ -(1:burn.in), ]
# store.x.2 <- store.x.2[ -(1:burn.in), ]

#-- -- --
#plotting

plot.ts(store.x.1[,1] - store.x.2[,1])

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

#-- -- --
#plotting the differences

library(ggplot2)
library(latex2exp)

len <- 1000

p1 <- ggplot() +
  geom_line(aes(x = 1:len, y = (store.x.1[,1] - store.x.2[,1])[1:len]),
            col = "darkslategray4", lwd = 0.75) +
  labs(x = "Iterations",
       y = TeX(r"($X_1 - X_2$ )",),
       col = "Index") +
  ggtitle("First component") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14))

p2 <- ggplot() +
  geom_line(aes(x = 1:len, y = (store.x.1[,2] - store.x.2[,2])[1:len]),
            col = "darkslategray4", lwd = 0.75) +
  labs(x = "Iterations",
       y = TeX(r"($X_1 - X_2$ )",),
       col = "Index") +
  ggtitle("Second component") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14))

p3 <- ggplot() +
  geom_line(aes(x = 1:len, y = (store.x.1[,3] - store.x.2[,3])[1:len]),
            col = "darkslategray4", lwd = 0.75) +
  labs(x = "Iterations",
       y = TeX(r"($X_1 - X_2$ )",),
       col = "Index") +
  ggtitle("Third component") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14))

p4 <- ggplot() +
  geom_line(aes(x = 1:len, y = (store.x.1[,4] - store.x.2[,4])[1:len]),
            col = "darkslategray4", lwd = 0.75) +
  labs(x = "Iterations",
       y = TeX(r"($X_1 - X_2$ )",),
       col = "Index") +
  ggtitle("Fourth component") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14))

p5 <- ggplot() +
  geom_line(aes(x = 1:len, y = (store.x.1[,5] - store.x.2[,5])[1:len]),
            col = "darkslategray4", lwd = 0.75) +
  labs(x = "Iterations",
       y = TeX(r"($X_1 - X_2$ )",),
       col = "Index") +
  ggtitle("Fifth component") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14))

p6 <- ggplot() +
  geom_line(aes(x = 1:len, y = (store.x.1[,6] - store.x.2[,6])[1:len]),
            col = "darkslategray4", lwd = 0.75) +
  labs(x = "Iterations",
       y = TeX(r"($X_1 - X_2$ )",),
       col = "Index") +
  ggtitle("Sixth component") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14))


library(patchwork)

p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)


