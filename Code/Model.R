
# --------- Packages -------------
library(deSolve)
library(extraDistr)
library(EnvStats)
library(truncdist)
library(tidyr)
library(ggplot2)

# --------- Function -------------
f_alp <- function(k, alp){
  alp[4] * (alp[3] / (1 + exp(2 * log(99) / alp[2] * (k - alp[1] - alp[2] / 2))) + 1 - alp[3])
}

eqn <- function(time, init, para, N){
  
  dS <- - para[1] * init[2] * init[1]  / N
  dE <- para[1] * init[2] * init[1] / N  + init[6] * para[3] - para[2] * init[2]
  dI <- para[4] * init[6] - init[3] * sum(para[5:6]) + para[9] * init[2]
  dR <- para[5] * init[3]
  dD <- para[6] * init[3]
  dPE <- - init[6] * (para[3] + para[4]) + (para[2] - para[9]) * init[2]
  
  return(list(c(dS, dE, dI, dR, dD, dPE)))
}

time_l <- function(k, Y, para_1, para, N, a, b, limit_P){
  para[1] <- para_1[k]
  out <- as.numeric(ode(y = Y[k,], times = 0:1, eqn, parms = para, N = N, atol = 1)[2,-1]) + 0.01
  Lik <- sum(dtpois(floor(Y[k+1,a:5]), out[a:5], -10^(-10), b = N, log = T))
  
  if(is.finite(Lik) == F){
    Lik <- -10^(10)
  }
  
  return(Lik)
}

l <- function(para, Y, E, P, N, a, b, para_1, limit_P, mark){
  Y[,2] <- E
  Y[,6] <- P
  Y[,1] <- N - rowSums(Y[,-1]) 
  
  logLik_Y <- sapply(mark, time_l, Y = Y, para_1 = para_1, para = para, N = N, a = a, b = b, limit_P = limit_P)
  return(sum(logLik_Y))
}

gibbs <- function(para, Y, N, time_change, dead, prior, E_max){ 
  
  E <- para[[2]]
  P <- para[[3]]
  alp <- para[[4]]
  q <- para[[5]]
  max_l <- para[[6]]
  para <- para[[1]]
  
  para_1 <- f_alp(1:(nrow(Y)-1), alp)
  
  s <- 5
  para[3] <- 0
  para[4] <- 0
  para[9] <- para[2]
  
  for(i in 5){
    para_t <- para
    para_t[i] <- rlnormTrunc(1, meanlog = log(para[i]), sdlog = log(1.1), min = 10^(-10), max = 1)
    b <- l(para_t, Y, E, P, N, 2, s, para_1, limit_P, 1:(nrow(Y)-1)) 
    r <- b - max_l 
    U <- log(runif(1))
    if(U < r){
      para <- para_t
      max_l <- b
    }
  }
  
  if(dead == T){
    for(i in 6){
      para_t <- para
      para_t[i] <- rlnormTrunc(1, meanlog = log(para[i]), sdlog = log(1.1), min = 10^(-10), max = 1)
      b <- l(para_t, Y, E, P, N, 2, s, para_1, limit_P, 1:(nrow(Y)-1)) 
      r <- b - max_l 
      U <- log(runif(1))
      if(U < r){
        para <- para_t
        max_l <- b
      }
    }
  }else{
    para[6] <- 0
  }
  
  Y[,6] <- P
  Test <- runif(1) < 0.5
  if(Test == T){
    Y[1,2] <- rtnorm(1, E[1], prior[3], 1, prior[2])
    Y[1,1] <- N - sum(Y[1,-1])
    
    mark <- 1:(nrow(Y)-1)
    
    for(k in mark){
      para[1] <- para_1[k]
      out <- as.numeric(ode(y = as.numeric(Y[k,]), times = 0:1, eqn, parms = para, N = N)[2,-1])
      Y[k + 1, 2] <- rtpois(1, out[2], a = 10^(-10) , b = E_max)
      Y[k + 1, 1] <- N - sum(Y[k + 1,-1])
    }
    
    E_t <- Y[,2]
    r <- l(para, Y, E_t, P, N, 3, s, para_1, limit_P, 1:(nrow(Y)-1)) - 
      l(para, Y, E, P, N, 3, s, para_1, limit_P, 1:(nrow(Y)-1))
    U <- log(runif(1))
    if(U < r){
      E <- E_t
    }
    
  }else{
    for(k in 1){
      E_t <- E
      E_t[k] <- rtnorm(1, E[1], prior[3], 1, prior[2])
      r <- l(para, Y, E_t, P, N, 2, s, para_1, limit_P, 1)-l(para, Y, E, P, N, 2, s, para_1, limit_P, 1)
      U <- log(runif(1))
      if(U < r){
        E <- E_t
      }
    }
    
    for(k in 2:(nrow(Y)-1)){
      E_t <- E
      E_t[k] <- rtnorm(1, E[k], prior[3],  10^(-10), E_max)
      r <- l(para, Y, E_t, P, N, 2, s, para_1, limit_P, (k-1):k)-l(para, Y, E, P, N, 2, s, para_1, limit_P, (k-1):k)
      U <- log(runif(1))
      if(U < r){
        E <- E_t
      }
    }
    
    for(k in nrow(Y)){
      E_t <- E
      E_t[k] <- rtnorm(1, E[k], prior[3],  10^(-10), E_max)
      r <- l(para, Y, E_t, P, N, 2, s, para_1, limit_P, (k-1))-l(para, Y, E, P, N, 2, s, para_1, limit_P, (k-1))
      U <- log(runif(1))
      if(U < r){
        E <- E_t
      }
    }
    
    max_l <- l(para, Y, E, P, N, 2, s, para_1, limit_P, 1:(nrow(Y)-1))
  }
  
  if(time_change == T){
    alp_t <- alp
    alp_t[1] <- rtnorm(1, alp[1], 1, 9, 50)
    b <- l(para, Y, E, P, N, 2, s, f_alp(1:(nrow(Y)-1), alp_t), limit_P, 1:(nrow(Y)-1))
    r <- b - max_l
    U <- log(runif(1))
    if(U < r){
      alp <- alp_t
      max_l <- b
    }
    
    alp_t <- alp
    alp_t[2] <- rtnorm(1, alp[2], 1, 1, 30)
    b <- l(para, Y, E, P, N, 2, s, f_alp(1:(nrow(Y)-1), alp_t), limit_P, 1:(nrow(Y)-1))
    r <- b - max_l
    U <- log(runif(1))
    if(U < r){
      alp <- alp_t
      max_l <- b
    }
  }
  
  alp_t <- alp
  alp_t[4] <- rlnormTrunc(1, meanlog = log(alp[4]), sdlog = log(1.1), min = 0.01)
  b <- l(para, Y, E, P, N, 2, s, f_alp(1:(nrow(Y)-1), alp_t), limit_P, 1:(nrow(Y)-1)) 
  r <- b - max_l 
  U <- log(runif(1))
  if(U < r){
    alp <- alp_t
    max_l <- b
  }
  
  return(list(para, E, P, alp, q, max_l))
}

para_est <- function(dat_tol, N, time_change){ 
  
  if(sum(dat_tol[,5]) == 0){
    dead <- F
  }else{
    dead <- T
  }
  
  E_max <- Inf
  E_prior <-  as.numeric(dat_tol[,3]) / 2
  
  dat_tol[,3] <- dat_tol[,3] - dat_tol[,4] - dat_tol[,5]
  dat_tol <- dat_tol + 0.001
  
  dead_prior <- mean((dat_tol[-1,5] - dat_tol[-nrow(dat_tol), 5]) / dat_tol[-nrow(dat_tol), 3]) + 0.01
  revo_prior <- mean((dat_tol[-1,4] - dat_tol[-nrow(dat_tol), 4]) / dat_tol[-nrow(dat_tol), 3]) + 0.01
  
  b <- rep(0.1, 9)
  b[5] <- revo_prior
  b[6] <- dead_prior
  b[2] <- 1 / 5
  b[9] <- 0
  
  if(time_change == T){
    a <- c(10, 2, 1, 2 * b[2])
  }else{
    a <- c(10, 1, 0, 0.5)
  }
  
  prior <- c(1, 4 * dat_tol[1,3], (dat_tol[1,3] + 1) * 0.2)
  para <- list()
  para[[1]] <- list(b,  E_prior, c(dat_tol[1,6] * 0.5, dat_tol[-1,6] * 0.5), a, c(0.2, 0.4), -Inf)
  
  k <- 10000
  for(i in 2:k){
    para[[i]] <- gibbs(para[[i-1]], dat_tol, N, time_change, dead, prior, E_max) 
    cat(c(i, signif(para[[i]][[1]][c(3:4)]), para[[i]][[2]][1:2], para[[i]][[3]][1:2], para[[i]][[4]][-3]), "\n")
  }
  
  mark <- 5001:10000
  
  para <- lapply(mark, function(k){
    para[[k]]
  }) 
  
  return(para)
}

pred <- function(para, dat_tol, d, N){ 
  
  dat_tol[,3] <- dat_tol[,3] - dat_tol[,4] - dat_tol[,5]
  
  pre_arr <- array(0, c(length(para), 6, d + nrow(dat_tol)))
  
  for(j in 1:length(para)){
    Y <- dat_tol
    Y[,2] <- para[[j]][[2]]
    Y[,6] <- para[[j]][[3]]
    Y[,1] <- N - rowSums(Y[,-1]) 
    
    pre <- Y[1,]
    
    for(k in 2:(nrow(dat_tol))){
      a <- para[[j]][[1]]
      a[1] <- f_alp(k-1, para[[j]][[4]])
      out <- as.numeric(ode(y = as.numeric(Y[k-1,]), times = 0:1, eqn, parms = a, N = N)[2,-1])
      pre <- cbind(pre, rep(0, 6))
      pre[,k] <- out
    }
    
    for(k in (nrow(dat_tol)+1):(nrow(dat_tol)+d)){
      a <- para[[j]][[1]]
      a[1] <- f_alp(k-1, para[[j]][[4]])
      out <- as.numeric(ode(y = as.numeric(Y[k-1,]), times = 0:1, eqn, parms = a, N = N)[2,-1])
      pre <- cbind(pre, rep(0, 6))
      pre[,k] <- out
      
      Y <- rbind(Y, rep(0, 6))
      Y[k,] <- out
    }
    
    pre_arr[j,,] <- pre
  }
  
  return(pre_arr)
}


TrainModel <- function(dat, N, d, time_change) {
  
  dat_tol <- dat
  para <- para_est(dat, N, time_change)
  pre_arr <- pred(para, dat_tol, d, N) 
  
  para_mean <- lapply(1:5, function(j){
    rowMeans(sapply(1:length(para), function(k){
      para[[k]][[j]]
    }))
  })
  
  E_Quan <- sapply(1:(d + nrow(dat_tol)), function(j){
    return(c(quantile(pre_arr[,2,j], c(0.025)), quantile(pre_arr[,2,j], c(0.5)), quantile(pre_arr[,2,j], c(0.975))))
  })
  
  I_Quan <- sapply(1:(d + nrow(dat_tol)), function(j){
    return(c(quantile(pre_arr[,3,j], c(0.025)), quantile(pre_arr[,3,j], c(0.5)), quantile(pre_arr[,3,j], c(0.975))))
  })
  
  R_Quan <- sapply(1:(d + nrow(dat_tol)), function(j){
    a <- pre_arr[,4,j]
    return(c(quantile(a, c(0.025)), quantile(a, c(0.5)), quantile(a, c(0.975))))
  })
  
  D_Quan <- sapply(1:(d + nrow(dat_tol)), function(j){
    a <- pre_arr[,5,j]
    return(c(quantile(a, c(0.025)), quantile(a, c(0.5)), quantile(a, c(0.975))))
  })
  
  P_Quan <- sapply(1:(d + nrow(dat_tol)), function(j){
    a <- pre_arr[,6,j] 
    return(c(quantile(a, c(0.025)), quantile(a, c(0.5)), quantile(a, c(0.975))))
  })
  
  All_Quan <- sapply(1:(d + nrow(dat_tol)), function(j){
    a <- pre_arr[,5,j] + pre_arr[,4,j] + pre_arr[,3,j]
    return(c(quantile(a, c(0.025)), quantile(a, c(0.5)), quantile(a, c(0.975))))
  })
  
  incu <- sapply(1:length(para), function(k){
    1 / para[[k]][[1]][2]
  })
  incu_Quan <- quantile(incu, c(0.025, 0.5, 0.975))
  
  cure <- sapply(1:length(para), function(k){
    1 / sum(para[[k]][[1]][c(5,6)])
  })
  cure_Quan <- quantile(cure, c(0.025, 0.5, 0.975))
  
  a <- sapply(1:length(para), function(k){
    para[[k]][[4]][c(4)] / para[[k]][[1]][c(2)]
  })
  R_0_Quan <- quantile(a, c(0.025, 0.5, 0.975))
  
  error_Quan <- NULL
  Dia_Quan <- NULL
  
  result <- list(para_mean = para_mean, 
                 E_Quan = E_Quan, I_Quan = I_Quan, 
                 D_Quan = D_Quan, R_Quan = R_Quan, 
                 P_Quan = P_Quan, All_Quan = All_Quan, 
                 incu_Quan = incu_Quan, Dia_Quan = Dia_Quan,
                 cure_Quan = cure_Quan, 
                 R_0_Quan = R_0_Quan, error_Quan = error_Quan,
                 para = para, pre_arr = pre_arr)
  
  return(list(dat = dat, N = N, result = result))
}

