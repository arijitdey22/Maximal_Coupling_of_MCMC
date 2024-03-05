# library(scales)
# library(mvtnorm)
# 
# ###=============================================================================
# ###Comparison between meeting time of Algorithm 2 and Algorithm 3
# 
# #Target: Multivariate Normal with mean = rep(2,d) and sigma:
#     #Case I: 2 * diag(d)
#     #Case II: CAR structure. We do this multiple times for varying rho.
# 
# rm(list = ls())
# 
# ##--------------------------------------
# #necessary functions for algorithm 3
# 
# fun.a.x.y <- function(x.1,y.1, d, target.sigma){
#   min(1, dmvnorm(y.1, mean = rep(2,d), sigma = target.sigma) / dmvnorm(x.1, mean = rep(2,d), sigma = target.sigma))
# }
# 
# fun.q.x.y <- function(x.2, y.2, cand.sd){
#   dmvnorm(y.2, x.2, cand.sd^2 * diag(d))
# }
# 
# fun.f.x.y <- function(x.3, y.3, cand.sd, d, target.sigma){
#   fun.a.x.y(x.3,y.3, d, target.sigma) * fun.q.x.y(x.3, y.3, cand.sd)
# }
# 
# ##--------------------------------------
# 
# d.vec <- c(2,3,4,5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30)
# 
# store.algo.2 <- numeric(length = length(d.vec))
# store.algo.3 <- numeric(length = length(d.vec))
# 
# for (i.d in 1:length(d.vec)){
#   
#   d <- d.vec[i.d]                                       #dimension 
#   
#   target.mean <- rep(2,d)                               #target distribution
#   target.sigma <- 2 * diag(d)
#   # foo <- matrix(rep(1:d, each = d), ncol = d)
#   # target.sigma <- 0.5 ^ abs(foo - t(foo))
#   
#   ##--------------------------------------
#   ##algorithm 02:
#   
#   cand.sd <- 1                                    #standard deviation of the proposal distribution
#   
#   reps <- 25                                      #number of replications for stabilizing results
#   count.iter.2.vec <- numeric(length = reps)      #storage for counts
#   
#   for (i.reps in 1:reps){        #loop for reps
#     
#     #set.seed(i.reps)
#     
#     x.cur.2 <- runif(d)            #initial values 
#     y.cur.2 <- runif(d)
#     
#     y.accept.alg.2 <- 0            #marker variable for acceptance
#     count.iter.2 <- 0                #count of the first meeting
#     
#     while (y.accept.alg.2 == 0){
#       
#       count.iter.2 <- count.iter.2 + 1
#       
#       #proposal of x
#       x.cand.2 <- rmvnorm(1, mean = x.cur.2, sigma = cand.sd^2 * diag(d))
#       U <- runif(1)
#       val.U <- dmvnorm(x.cand.2, mean = y.cur.2, sigma = cand.sd^2 * diag(d)) / dmvnorm(x.cand.2, mean = x.cur.2, sigma = cand.sd^2 * diag(d))
#       
#       #proposal of y
#       if (U <= val.U){
#         y.cand.2 <- x.cand.2
#       }else{
#         accept.2 = 0
#         while(accept.2 == 0){
#           y.curl.2 <- rmvnorm(1, mean = y.cur.2, sigma = cand.sd^2 * diag(d))
#           V <- runif(1)
#           V.val <- dmvnorm(y.curl.2, mean = x.cur.2, sigma = cand.sd^2 * diag(d)) / dmvnorm(y.curl.2, mean = y.cur.2, sigma = cand.sd^2 * diag(d))
#           
#           if(V > V.val){
#             y.cand.2 <- y.curl.2
#             accept.2 = 1
#           }
#         }
#       }
#       
#       #acceptance probabilities
#       a.xx <- fun.a.x.y(x.cur.2, x.cand.2, d, target.sigma)
#       a.yy <- fun.a.x.y(y.cur.2, y.cand.2, d, target.sigma)
#       
#       B.x <- rbinom(1, 1, min(1, a.xx))
#       B.y <- rbinom(1, 1, min(1, a.yy))
#       
#       if ( ( sum(x.cand.2 != y.cand.2) == 0 ) & (B.x == B.y) ){
#         y.accept.alg.2 <- 1
#       }
#       
#     }
#     
#     #aesthetics
#     print(paste0("Algorithm 2; d = ", d.vec[i.d], ", rep = ", i.reps))
#     
#     count.iter.2.vec[i.reps] <- count.iter.2
#   }
#   
#   ##--------------------------------------
#   ##algorithm 03:
#   
#   cand.sd <- 1                                    #standard deviation of the proposal distribution
# 
#   reps <- 20                                      #number of replications for stabilizing results
#   count.iter.3.vec <- numeric(length = reps)      #storage for counts
# 
#   for (i.reps in 1:reps){        #loop for reps
# 
#     #set.seed(i.reps)
# 
#     x.cur.3 <- runif(d)            #initial values
#     y.cur.3 <- runif(d)
# 
#     y.accept.alg.3 <- 0            #marker variable for acceptance
#     count.iter.3 <- 0                #count of the first meeting
# 
#     while (y.accept.alg.3 == 0){
# 
#       count.iter.3 <- count.iter.3 + 1
# 
#       #proposal of x
#       x.cand.3 <- rmvnorm(1, mean = x.cur.3, sigma = cand.sd^2 * diag(d))
# 
#       #accepting x
#       a.xx.3 <- fun.a.x.y(x.cur.3, x.cand.3, d, target.sigma)
#       B.x.3 <- rbinom(1, 1, min(1, a.xx.3))
#       x.prev.3 <- x.cur.3
#       x.cur.3 <- B.x.3 * x.cand.3 + (1-B.x.3) * x.cur.3
# 
#       #accepting y
#       f.x.X.3 <- fun.f.x.y(x.prev.3, x.cur.3, cand.sd, d, target.sigma)
#       f.y.X.3 <- fun.f.x.y(y.cur.3, x.cur.3, cand.sd, d, target.sigma)
#       U <- runif(1)
# 
#       if ( (sum(x.cur.3 == x.prev.3) == 0) & (U <=  (f.x.X.3 / f.y.X.3))){
#         y.cur.3 <- x.cur.3
#         y.accept.alg.3 <- 1
#       }else{
# 
#         accpet.y.3 <- 0
# 
#         while (accpet.y.3 == 0){
# 
#           y.curl.3 <- rmvnorm(1, mean = y.cur.3, sigma = cand.sd^2 * diag(d))
# 
#           if (sum(y.curl.3 != y.cur.3) == 0){
#             y.cur.3 = y.curl.3
#             accpet.y.3 = 1
# 
#           }else{
#             f.y.y.curl.3 <- fun.f.x.y(y.cur.3, y.curl.3, cand.sd, d, target.sigma)
#             f.x.y.curl.3 <- fun.f.x.y(x.cur.3, y.curl.3, cand.sd, d, target.sigma)
#             V <- runif(1)
# 
#             if (V < (f.x.y.curl.3 / f.y.y.curl.3)){
#               y.cur.3 = y.curl.3
#               accpet.y.3 = 1
# 
#             }
#           }
# 
#         }
# 
#       }
# 
#     }
# 
#     #aesthetics
#     print(paste0("Algorithm 3; d = ", d.vec[i.d], ", rep = ", i.reps))
# 
#     count.iter.3.vec[i.reps] <- count.iter.3
#   }
#   
#   store.algo.2[i.d] <- mean(count.iter.2.vec)
#   store.algo.3[i.d] <- mean(count.iter.3.vec)
#   
# }
# 
# ##--------------------------------------
# 
# save(store.algo.2, store.algo.3, file = "03_Comparison.1.Rdata")

load("03_Comparison.1.Rdata")

library(ggplot2)

df.1 <- data.frame(rep(d.vec,2), c(store.algo.2, store.algo.3),
                   c(rep("Algo 2", length(d.vec)), rep("Algo 3", length(d.vec)) ))
colnames(df.1) <- c("d", "n", "type")

ggplot(data = df.1) +
  geom_line(aes(x = d, y = n, col = type), lwd = 0.8) +
  labs(x = "Dimensions",
       y = "Meeting time",
       col = "") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 13),
        legend.position = "bottom")

###=============================================================================
###plot of d vs Gibbs maximal coupling time



