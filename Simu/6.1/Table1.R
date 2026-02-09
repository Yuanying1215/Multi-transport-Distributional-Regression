# ==============================================================================
# Script: Table1.R
# Description: Simulation study for the Single Predictor setting (Section 6.1).
#              This script evaluates the estimation performance of the MTDR model.
# Output: Reproduces the Mean and SD of estimation errors reported in Table 1.
#
# Method: Multi-transport Distributional Regression (MTDR)
# Key Algorithm: Iterative Back-fitting (Algorithm 1 in the paper)
# ==============================================================================

# --- Core Algorithm Functions ---
#' gpava_new():
#' Generalized Pool Adjacent Violators Algorithm (gPAVA)
#' Solves the weighted isotonic regression sub-problem for updating transport maps.

#' Main_new():
#' Main Iterative Optimization (Algorithm 1)
#' Alternates between updating transport maps (via gPAVA) and weights (via optimization).
#' 
#' parameters:
#' @param a_true,m,n 
#' @param al,ar Lower truncation threshold/Upper truncation threshold
#' @param zero Small constant to prevent division by zero

gpava_new <- function (z,y,weights,s0,s1,solver = weighted.mean) 
{
  rule <- 1
  x <- z
  w <- weights
  merger <- c 
  if ((is.matrix(y)) || (is.data.frame(y))) 
    y <- c(as.data.frame(t(y)))
  if ((is.matrix(w)) || (is.data.frame(w))) 
    w <- c(as.data.frame(t(w)))
  w1 <- w
  y1 <- y
  n <- length(y)
  inds <- as.list(seq_len(n))
  
  vals <- mapply(function(y, w) {
    if (sum(w, na.rm = TRUE) > 0) {
      weighted.mean(y, w, na.rm = TRUE)
    } else {
      mean(y, na.rm = TRUE)  #Return the simple (arithmetic) mean when the sum of weights is zero.
    }
  }, y, w)
  if (rule == 2){vals <- pmin(pmax(vals,s0),s1)}
  combine <- function(i) {
    j <- i + 1L
    y[[i]] <<- merger(y[[i]], y[[j]])
    w[[i]] <<- c(w[[i]], w[[j]])
    
    vals[i] <<- solver(y[[i]], w[[i]])
    
    inds[[i]] <<- c(inds[[i]], inds[[j]])
    keep <- seq_len(n)[-j]
    y <<- y[keep]
    w <<- w[keep]
    vals <<- vals[keep]
    inds <<- inds[keep]
    n <<- n - 1L
  }
  i <- 1L
  repeat {
    if (i < n) {
      if ((vals[i] > vals[i + 1])) {
        combine(i)
        while ((i > 1L) && (vals[i - 1L] > vals[i])) {
          combine(i - 1L)
          i <- i - 1L
        }
      }
      else i <- i + 1L
    }
    else break
  }
  if (rule == 1){vals <- pmin(pmax(vals,s0),s1)}
  if (rule != 0){
    vals[1] <- s0
    vals[n] <- s1
  }
  yfit <- rep.int(vals, sapply(inds, length))
  result <- list(x = yfit, z = x, y = y1, w = w1, solver = solver, 
                 call = match.call())
  class(result) <- "gpava"
  result
}
# lossfun <- function(a0,T0,T1,T01,T0nu,x,w0,N){
#   (1/N)*sum((a0*T0(x) + (1-a0)*T1(T01) - T0nu)^2 * w0)
# }
# lossfun <- function(a0,T0,T1,T10,T1nu,x,w1,N){
#   (1/N)*sum((a0*T0(T10) + (1-a0)*T1(x) - T1nu)^2 * w1)
# }
lossfun <- function(a0,T0,T1,T01,T0nu,w0,T10,T1nu,x,w1,N){
  (1/N)*sum((a0*T0(T10) + (1-a0)*T1(x) - T1nu)^2 * w1)+
    (1/N)*sum((a0*T0(x) + (1-a0)*T1(T01) - T0nu)^2 * w0)
}
Main_new <- function(T0,T1,a0,obj_0,T10,T1nu,w1,T01,T0nu,w0,x,tol,
                     al,ar,zero,a_type,max_iter=500)
{
  m <- length(x)
  s0 <- min(x);s1 <- max(x)
  N <- ncol(T10)
  for (iter in 1:max_iter) {
    #T1
    y <- (T1nu - a0*T0(T10))/(1-a0)
    T1_grid <- gpava_new(x,y,weights = w1,s0,s1)$x 
    T1_new <- approxfun(x, T1_grid, method = "linear")
    #T0
    y <- (T0nu - (1-a0)*T1_new(T01))/a0
    T0_grid <- gpava_new(x,y,weights = w0,s0,s1)$x 
    T0_new <- approxfun(x, T0_grid, method = "linear")
    #alpha
    if (a_type == 0){
      a <- sum((T0nu - T1_new(T01))*(T0_new(x)-T1_new(T01))*w0)/
        sum(((matrix((T0_new(x)- T1_new(T01)),nrow=m))^2)*w0)
      if (a < al){a <- zero}
      if (a > ar){a <- 1 - zero}
    } else {a <- optimize(lossfun,c(0,1),T0_new,T1_new,T01,T0nu,w0,T10,T1nu,x,w1,N)$minimum
    }
    #objective_loss
    obj1 <- (1/N)*sum((a*T0_new(x) + (1-a)*T1_new(T01) - T0nu)^2 * w0)
    obj2 <- (1/N)*sum((a*T0_new(T10) + (1-a)*T1_new(x) - T1nu)^2 * w1)
    obj <- obj1 + obj2
    cat(sprintf("Iteration %d: Objective = %.12f\n", iter, obj))
    if (iter > 1 && abs(obj_0 - obj) < tol) {
      cat("Converged.\n")
      break
    }
    obj_0 <- obj;T1 <- T1_new;T0 <- T0_new;a0 <- a
    iter <- iter + 1
  }
  return(list(T0_res = T0_new,T1_res = T1_new,a_res = a,loss_res=obj/2))
}



library(parallel)
#library(quadprog)
library(pracma)

a_true <- 0.5 #0#0.25#0.75#1
b <- 1 - a_true
al <- 0.0001; ar <- 1-al #
zero <- 1e-6 #
zeta <- function(x, k) {
  if (k == 0) {
    return(x)
  }
  x - sin(pi * k * x) / (pi * abs(k))
}
T0_true <- function(x){
  zeta(x,4)
}
T1_true <- function(x){
  zeta(x,3)
}
nqSup <- 101
qSup <- seq(0, 1, length.out = nqSup)
t <- 101
x <- seq(0, 1, length.out = t)
configs <- list(
  c(50,200),
  c(50,400),
  c(200,200),
  c(200,400),
  c(400,200),
  c(400,400)
)
for (config in configs) {
  n <- config[1]
  m <- config[2]
  test <- 0.3*n
  for (a_type in c(1)){
    B <-  100
    seeds <- seq(1,B,1)
    results <- list()
    cl <- makeCluster(20,type="FORK")  
    results <- parLapply(cl,1:B,function(l){
      set.seed(seeds[l])
      alpha1 <- runif(n+test, 1, 5)
      beta1 <- runif(n+test, 1, 5)
      
      X1 <- sapply( seq_len(n+test), function (i) {
        qbeta(runif(m), alpha1[i],beta1[i])
      })
      X1q <- sapply(seq_len(n+test), function (i) {
        quantile(X1[,i],qSup)
      })
      ecdf_X1 <- apply(X1q,2,function(col){
        approxfun(col, qSup, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      eqf_X1 <- apply(X1q,2,function(col){
        approxfun(qSup, col, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      X0 <- sapply(seq_len(n+test), function (i) {
        qunif(runif(m))
      })
      X0q <- sapply(seq_len(n+test), function (i) {
        quantile(X0[,i],qSup)
      })
      ecdf_bar <- apply(X0q,2,function(col){
        approxfun(col, qSup, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      eqf_bar <- apply(X0q,2,function(col){
        approxfun(qSup, col, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      # ecdf_bar <- lapply(1:(n+test), function(i) {
      #     function(x){x}
      # })
      # eqf_bar <- lapply(1:(n+test), function(i) {
      #     function(x){x}
      # })
      
      u <- sapply(seq_len(n+test), function (i) {
        runif(m)
      })
      X0r = sapply(seq_len(n+test), function (i) {
        qunif(u[,i])
      })
      X1r = sapply(seq_len(n+test), function (i) {
        qbeta(u[,i], alpha1[i],beta1[i])
      })
      Y_res <- b*T0_true(X0r) + a_true*T1_true(X1r)
      values <- setdiff(-20:20, c(0))
      K <- sample(values,n,replace = T)
      Y_sample <- Y_res
      # Y_sample[,1:n] <- sapply(seq_len(n), function(i){
      #     zeta(Y_sample[,i],K[i])
      # })
      for (i in 1:n) {
        Y_sample[, i] <- zeta(Y_sample[, i], K[i])
      }
      Y <- sapply(seq_len(n+test), function (i) {
        quantile(Y_sample[,i],qSup)
      })
      ecdf_Y <- apply(Y[,1:n],2,function(col){
        approxfun(col, qSup, method = "linear", 
                  yleft = 0, yright = 1, ties = "ordered")
      })
      eqf_Y <- apply(Y[,1:n],2,function(col){
        approxfun(qSup, col, method = "linear", 
                  yleft = 0, yright = 1, ties = "ordered")
      })
      T10 <- sapply(seq_len(n+test), function(i){
        eqf_bar[[i]](ecdf_X1[[i]](x))
      })
      T1nu <- sapply(seq_len(n), function(i){
        eqf_Y[[i]](ecdf_X1[[i]](x))
      })
      w1 <- sapply(seq_len(n+test), function(i){
        ecdf_X1[[i]](x)
      })
      w1 <- w1 - rbind(0,w1[-nrow(w1),])
      T01 <- sapply(seq_len(n+test), function(i){
        eqf_X1[[i]](ecdf_bar[[i]](x))
      })
      T0nu <- sapply(seq_len(n), function(i){
        eqf_Y[[i]](ecdf_bar[[i]](x))
      })
      w0 <- sapply(seq_len(n+test), function(i){
        ecdf_bar[[i]](x)
      })
      w0 <- w0 - rbind(0,w0[-nrow(w0),])
      T1 <- approxfun(x,x,method="linear")
      T0 <- approxfun(x,x,method="linear")
      
      #####
      if (a_type==0){
        a0 <- sum((T0nu - T1(T01))*(T0(x)-T1(T01))*w0)/
          sum(((matrix((T0(x)- T1(T01)),nrow=t))^2)*w0)
        if (a0 < al){a0 <- zero}
        if (a0 > ar){a0 <- 1 - zero}
      }else{
        a0 <- optimize(lossfun,c(0,1),T0=T0,T1=T1,T01=T01[,1:n],T0nu=T0nu,w0=w0[,1:n],T10=T10[,1:n],T1nu=T1nu,x=x,w1=w1[,1:n],N=n)$minimum
      }
      
      obj1_0 <- (1/n)*sum((a0*T0(x) + (1-a0)*T1(T01[,1:n]) - T0nu)^2 * w0[,1:n])
      obj2_0 <- (1/n)*sum((a0*T0(T10[,1:n]) + (1-a0)*T1(x) - T1nu)^2 * w1[,1:n])
      obj_0 <- obj1_0 + obj2_0
      cat(sprintf("Objective 0: Objective = %.12f\n", obj_0))
      
      result <- Main_new(T0,T1,a0,obj_0,T10[,1:n],T1nu,w1[,1:n],T01[,1:n],T0nu,w0[,1:n],x,tol=1e-9,
                         al,ar,zero,a_type)
      a_res <- 1 - result$a_res
      T0_res <- result$T0_res
      T1_res <- result$T1_res
      
      # RMSE_test0 <- sqrt((1/test)*sum(((1-a_res)*T0_res(x)
      #                      + (a_res)*T1_res(T01[,(n+1):(n+test)]) 
      #                      - b*T0_true(x)-a_true*T1_true(T01[,(n+1):(n+test)]))^2*w0[,(n+1):(n+test)]))
      # RMSE_test1 <- sqrt((1/test)*sum(((1-a_res)*T0_res(T10[,(n+1):(n+test)])
      #                      + (a_res)*T1_res(x) 
      #                      - b*T0_true(T10[,(n+1):(n+test)])-a_true*T1_true(x))^2*w1[,(n+1):(n+test)]))
      RMSE_test <- sqrt((1/test)*sum((((1-a_res)*T0_res(T10[,(n+1):(n+test)])
                                       + (a_res)*T1_res(x) 
                                       - b*T0_true(T10[,(n+1):(n+test)])-a_true*T1_true(x))^2*w1[,(n+1):(n+test)]+
                                        ((1-a_res)*T0_res(x)
                                         + (a_res)*T1_res(T01[,(n+1):(n+test)]) 
                                         - b*T0_true(x)-a_true*T1_true(T01[,(n+1):(n+test)]))^2*w0[,(n+1):(n+test)])/2))
      results[[l]] <- c(abs(a_res-a_true),sqrt(trapz(x, (T0_res(x) - T0_true(x))^2)),
                        sqrt(trapz(x, (T1_res(x) - T1_true(x))^2)),RMSE_test)
      
    })
    stopCluster(cl)
    
    a0_res =  mean(matrix(unlist(results),nrow = 4)[1,])
    a0_sd = sd(matrix(unlist(results),nrow = 4)[1,])
    L2_T0 = mean(matrix(unlist(results),nrow = 4)[2,])
    T0_sd = sd(matrix(unlist(results),nrow = 4)[2,])
    L2_T1 = mean(matrix(unlist(results),nrow = 4)[3,])
    T1_sd = sd(matrix(unlist(results),nrow = 4)[3,])
    AWD = mean(matrix(unlist(results),nrow = 4)[4,])
    AWD_sd = sd(matrix(unlist(results),nrow = 4)[4,])
    
    cat(sprintf("RMSE_test,a0,L2_T0,L2_T1: %.3f(%.3f) %.3f(%.3f) %.3f(%.3f) %.3f(%.3f)\n", 
                AWD,AWD_sd,a0_res,a0_sd,L2_T0,T0_sd,L2_T1,T1_sd))
    
    
  }
}



