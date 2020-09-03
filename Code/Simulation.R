# --------------- Packages ----------------
library(nCov2019)
library(deSolve)
library(extraDistr)
library(EnvStats)
library(truncdist)
library(tidyr)
library(ggplot2)

# --------------- Functions ----------------

pred <- function(para, dat_tol, d, N){

  dat_tol[,3] <- dat_tol[,3] - dat_tol[,4] - dat_tol[,5]

  pre_arr <- array(0, c(length(para), 6, d + 5))

  for(j in 1:length(para)){
    Y <- dat_tol
    Y[,2] <- para[[j]][[2]]
    Y[,6] <- para[[j]][[3]]
    Y[,1] <- N - rowSums(Y[,-1])
    pre_out <- matrix(0, nrow = d + 5, ncol = 6)
    pre_out[1:6, ] <- Y[1:6, ]

    for(k in 7:d){
      a <- para[[j]][[1]]
      a[1] <- f_alp(k-1, para[[j]][[4]])
      out <- as.numeric(ode(y = as.numeric(pre_out[k-1, ]), times = (k-1):k, eqn, parms = a, N = N)[2,-1])
      pre_out[k,] <- out
    }

    pre_arr[j,,] <- t(pre_out)
  }

  return(pre_arr)
}

change_d <- function(para, d) {
  for(i in 1:length(para)) {
    para[[i]][[4]][1] <- para[[i]][[4]][1] + d
  }
  return(para)
}

my_fun <- function(pre_arr, d) {
  E_Quan <- sapply(1:(d), function(j){
    quantile(pre_arr[,2,j], c(0.025, 0.5, 0.975))
  })

  I_Quan <- sapply(1:(d), function(j){
    quantile(pre_arr[,3,j], c(0.025, 0.5, 0.975))
  })

  All_Quan <- sapply(1:(d), function(j){
    quantile(pre_arr[,3,j] + pre_arr[,4,j] + pre_arr[,5,j], c(0.025, 0.5, 0.975))
  })

  R_Quan <- sapply(1:(d), function(j){
    quantile(pre_arr[,4,j], c(0.025, 0.5, 0.975))
  })

  D_Quan <- sapply(1:(d), function(j){
    quantile(pre_arr[,5,j], c(0.025, 0.5, 0.975))
  })

  P_Quan <- sapply(1:(d), function(j){
    quantile(pre_arr[,6,j], c(0.025, 0.5, 0.975))
  })

  result <- list(E_Quan = E_Quan, I_Quan = I_Quan,
                 D_Quan = D_Quan, R_Quan = R_Quan,
                 P_Quan = P_Quan, All_Quan = All_Quan)

  return(result)
}

formatter <- function(x) {
  if(x < 1e5) {
    return(x)
  }
  level <- floor(log10(x))
  temp <- round(x / 10^level, digits = 2)
  return(paste0(temp, " %*% 10^", level))
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


# ------------------- Simulation ----------------------
# start at 1.19
load("Data/result_sz.rda")
source("Code/Model.R")

(dat_tol <- train_result$dat)
para <- train_result$result$para

f_alp <- function(k, alp){
  alp[4] * (alp[3] / (1 + exp(2 * log(99) / alp[2] * (k - alp[1] - alp[2] / 2))) + 1 - alp[3])
}

dat1 <- load_nCov2019()
dat_city <- dat1$data
dat_city_sz <- dat_city[dat_city$city == "Shenzhen", ]
dat_city_sz <- na.omit(dat_city_sz)

ind = 42
plot(1:ind, dat_city_sz[1:ind, "cum_confirm"], ylim = c(0, 1000))

train_result$result$para_mean[[4]][1]
para_mean <- train_result$result$para_mean

pre_arr <- pred(para0, dat_tol = train_result$dat, 100, N = train_result$N)
result2 <- my_fun(pre_arr, ind)
result_output <- result2$All_Quan[, ]

# Plot
plot(1:ind, dat_city_sz[1:ind, "cum_confirm"], ylim = c(0, 3000))
lines(1:ind, result2$All_Quan[2, ][1:ind])

d = 100
for (d_late in 1:5) {
  pre_arr <- pred(change_d(para0, d_late), dat_tol = train_result$dat, 100, N = 13020000)
  result2 <- my_fun(pre_arr, ind)
  lines(1:ind, result2$All_Quan[2, 1:ind], col = d_late)
}
