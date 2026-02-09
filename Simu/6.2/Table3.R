# ==============================================================================
# Script: Table3.R
# Description: Simulation study for Multiple Predictors (Section 6.2).
#              Evaluates estimation performance with p=2 predictors.
# Output: Reproduces the parameter estimation errors reported in Table 3.
#
# Scenario: Two predictor distributions with different parameters.
# ==============================================================================

# --- Core Algorithm Functions ---
#' gpava_new():
#' Generalized Pool Adjacent Violators Algorithm (gPAVA)
#' Solves the weighted isotonic regression sub-problem for updating transport maps.

#' Main2():
#' Main Iterative Optimization (Algorithm 1)
#' Alternates between updating transport maps (via gPAVA) and weights (via optimization).
#' 
#' parameters:
#' @param a_true
#' @param N Number of distributions/subjects. 
#'         (Note: This is denoted as 'n' in the manuscript/paper)
#' @param n Number of samples per distribution.
#'         (Note: This is denoted as 'm' in the manuscript/paper) 
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
alpha_l <- function(T_list, OT_list, x, w0) 
{
  k <- length(T_list)       
  m <- length(x)            
  A <- matrix(0, nrow = k, ncol = k)
  b <- matrix(0, nrow = k, ncol = 1)
  for (i in 1:k) {
    Ti_x <- if (i == 1) T_list[[i]](x) else T_list[[i]](OT_list[[i-1]])
    b[i, 1] <- sum(rowMeans(Ti_x * OT_list[[k]]) * w0[, 1])
    for (j in i:k) {
      Tj_x <- if (j == 1) T_list[[j]](x) else T_list[[j]](OT_list[[j-1]])
      A[i, j] <- if (i==1 & j==1) sum((Ti_x * Tj_x) * w0[, 1]) else sum(rowMeans(matrix(Ti_x * Tj_x,nrow=m)) * w0[, 1])
      A[j, i] <- A[i, j]
    }
  }
  A[A == 0] <- 1e-4
  A_inv <- solve(A)
  test <- matrix((sum(A_inv %*% b) - 1) / sum(A_inv), nrow = k, ncol = 1)
  res <- A_inv %*% (b - test)
  return(res)
}
adjust_vector <- function(vec,al,ar,zero) 
{
  vec[vec < al] <- zero
  vec[vec > ar] <- 1-zero
  free_idx <- which(vec > 0 & vec < 1)
  fixed_sum <- sum(vec)
  vec <- vec/fixed_sum
  return(vec)
}
lossfun_multi <- function(a,T0,T1,T2,T01,T02,T0nu,T10,T12,T1nu,T20,T21,T2nu,x,w0,w1,w2,N){
  (1/N)*sum((a[1]*T0(x) + a[2]*T1(T01) + (1-a[1]-a[2])*T2(T02) - T0nu)^2 * w0)+
    (1/N)*sum((a[1]*T0(T10) + a[2]*T1(x) + (1-a[1]-a[2])*T2(T12) - T1nu)^2 * w1)+
    (1/N)*sum((a[1]*T0(T20) + a[2]*T1(T21) + (1-a[1]-a[2])*T2(x) - T2nu)^2 * w2)
}
Main2 <- function(T0,T1,T2,a,obj_0,T10,T12,T1nu,w1,T20,T21,T2nu,w2,
                  T01,T02,T0nu,w0,x,tol,al,ar,zero,a_type,max_iter=500)
{
  m <- length(x)
  s0 <- min(x);s1 <- max(x)
  N <- ncol(T10)
  for (iter in 1:max_iter) {
    #T0
    y <- (T0nu - a[2]*T1(T01) - a[3]*T2(T02))/a[1]
    T0_grid <- gpava_new(x,y,weights = w0,s0,s1)$x 
    T0_new <- approxfun(x, T0_grid, method = "linear")
    #T1
    y <- (T1nu - a[1]*T0_new(T10) - a[3]*T2(T12))/a[2]
    T1_grid <- gpava_new(x,y,weights = w1,s0,s1)$x 
    T1_new <- approxfun(x, T1_grid, method = "linear")
    #T2
    y <- (T2nu - a[1]*T0_new(T20) - a[2]*T1_new(T21))/a[3]
    T2_grid <- gpava_new(x,y,weights = w2,s0,s1)$x 
    T2_new <- approxfun(x, T2_grid, method = "linear")
    if (a_type == 0){
      anew <- alpha_l(list(T0_new,T1_new,T2_new),list(T01,T02,T0nu),x,w0)
    }else{
      aop <- optim(par=c(0,0),fn=lossfun_multi,T0=T0_new,T1=T1_new,
                   T2=T2_new,T01=T01,T02=T02,T0nu=T0nu,T10=T10,T12=T12,
                   T1nu=T1nu,T20=T20,T21=T21,T2nu=T2nu,x=x,w0=w0,w1=w1,
                   w2=w2,N=N)$par
      anew <- c(aop[1],aop[2],1-aop[1]-aop[2])
    }
    anew <- adjust_vector(anew,al,ar,zero)
    obj1 <- (1/N)*sum((anew[1]*T0_new(x) + anew[2]*T1_new(T01) + anew[3]*T2_new(T02) - T0nu)^2 * w0)
    obj2 <- (1/N)*sum((anew[1]*T0_new(T10) + anew[2]*T1_new(x) + anew[3]*T2_new(T12) - T1nu)^2 * w1)
    obj3 <- (1/N)*sum((anew[1]*T0_new(T20) + anew[2]*T1_new(T21) + anew[3]*T2_new(x) - T2nu)^2 * w2)
    obj <- obj1 + obj2 + obj3
    cat(sprintf("Iteration %d: Objective = %.12f\n", iter, obj))
    if (iter > 1 && abs(obj_0 - obj) < tol) {
      cat("Converged.\n")
      break
    }
    obj_0 <- obj;T0 <- T0_new;T1 <- T1_new;T2 <- T2_new;a <- anew
    iter <- iter + 1
  }
  return(list(T_res=list(T0_res = T0_new,T1_res = T1_new,T2_res = T2_new),
              a_res = anew,loss_res=obj/3))
}

library(parallel)
library(pracma)

a_true <- c(0.3,0.35,0.35)#c(0.2,0.2,0.6)#c(0,0.5,0.5)#c(0,1,0)
al <- 0; ar <- 1
zero <- 1e-6
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
T2_true <- function(x){
  zeta(x,-5)
}
nqSup <- 101
qSup <- seq(0, 1, length.out = nqSup)
m <- 101
x <- seq(0, 1, length.out = m)
configs <- list(
  c(50,200),
  c(50,400),
  c(200,200),
  c(200,400),
  c(400,200),
  c(400,400)
)
for (a_type in c(1)) {
  for (config in configs) {
    N <- config[1]
    n <- config[2]
    test <- 0.3*N
    B <-  100
    seeds <- seq(1,B,1)
    results <- list()
    cl <- makeCluster(20,type="FORK") 
    results <- parLapply(cl,1:B,function(l){
      set.seed(seeds[l])
      alpha1 <- runif(N+test, 1, 5)
      beta1 <- runif(N+test, 1, 5)
      alpha2 <- runif(N+test, 2,6)
      beta2 <- runif(N+test, 2,6)
      
      X1 <- sapply( seq_len(N+test), function (i) {
        qbeta(runif(n), alpha1[i],beta1[i])
      })
      X1q <- sapply(seq_len(N+test), function (i) {
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
      X2 <- sapply( seq_len(N+test), function (i) {
        qbeta(runif(n), alpha2[i],beta2[i])
      })
      X2q <- sapply(seq_len(N+test), function (i) {
        quantile(X2[,i],qSup)
      })
      ecdf_X2 <- apply(X2q,2,function(col){
        approxfun(col, qSup, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      eqf_X2 <- apply(X2q,2,function(col){
        approxfun(qSup, col, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      X0 <- sapply( seq_len(N+test), function (i) {
        qunif(runif(n))
      })
      X0q <- sapply(seq_len(N+test), function (i) {
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
      
      # ecdf_bar <- lapply(1:(N+test),function(i){
      #     function(x){x}
      # })
      # eqf_bar <- lapply(1:(N+test),function(i){
      #     function(x){x}
      # })
      #
      u <- sapply(seq_len(N+test), function (i) {
        runif(n)
      })
      X0r = sapply(seq_len(N+test), function (i) {
        qunif(u[,i])
      })
      X1r = sapply(seq_len(N+test), function (i) {
        qbeta(u[,i], alpha1[i],beta1[i])
      })
      X2r = sapply(seq_len(N+test), function (i) {
        qbeta(u[,i], alpha2[i],beta2[i])
      })
      Y_res <- a_true[1]*T0_true(X0r) + a_true[2]*T1_true(X1r) +a_true[3]*T2_true(X2r)
      #K <- sample(-100:100,N+test,replace = T) #K
      values <- setdiff(-20:20, c(0))
      K <- sample(values,N+test,replace = T)
      Y_sample <- sapply(seq_len(N+test), function(i){
        zeta(Y_res[,i],K[i])
      })
      Y <- sapply(seq_len(N+test), function (i) {
        quantile(Y_sample[,i],qSup)
      })
      ecdf_Y <- apply(Y[,1:N],2,function(col){
        approxfun(col, qSup, method = "linear", 
                  yleft = 0, yright = 1, ties = "ordered")
      })
      eqf_Y <- apply(Y[,1:N],2,function(col){
        approxfun(qSup, col, method = "linear", 
                  yleft = 0, yright = 1, ties = "ordered")
      })
      
      T10 <- sapply(seq_len(N+test), function(i){
        eqf_bar[[i]](ecdf_X1[[i]](x))
      })
      T12 <- sapply(seq_len(N+test), function(i){
        eqf_X2[[i]](ecdf_X1[[i]](x))
      })
      T1nu <- sapply(seq_len(N), function(i){
        eqf_Y[[i]](ecdf_X1[[i]](x))
      })
      w1 <- sapply(seq_len(N+test), function(i){
        ecdf_X1[[i]](x)
      })
      w1 <- w1 - rbind(0,w1[-nrow(w1),])
      T20 <- sapply(seq_len(N+test), function(i){
        eqf_bar[[i]](ecdf_X2[[i]](x))
      })
      T21 <- sapply(seq_len(N+test), function(i){
        eqf_X1[[i]](ecdf_X2[[i]](x))
      })
      T2nu <- sapply(seq_len(N), function(i){
        eqf_Y[[i]](ecdf_X2[[i]](x))
      })
      w2 <- sapply(seq_len(N+test), function(i){
        ecdf_X2[[i]](x)
      })
      w2 <- w2 - rbind(0,w2[-nrow(w2),])
      T01 <- sapply(seq_len(N+test), function(i){
        eqf_X1[[i]](ecdf_bar[[i]](x))
      })
      T02 <- sapply(seq_len(N+test), function(i){
        eqf_X2[[i]](ecdf_bar[[i]](x))
      })
      T0nu <- sapply(seq_len(N), function(i){
        eqf_Y[[i]](ecdf_bar[[i]](x))
      })
      w0 <- sapply(seq_len(N+test), function(i){
        ecdf_bar[[i]](x)
      })
      w0 <- w0 - rbind(0,w0[-nrow(w0),])
      
      T2 <- approxfun(x,x,method="linear")
      T1 <- approxfun(x,x,method="linear")
      T0 <- approxfun(x,x,method="linear")
      if (a_type==0){
        a <- alpha_l(list(T0,T1,T2),list(T01,T02,T0nu),x,w0)
      }else{
        aop <- optim(par=c(0,0),fn=lossfun_multi,T0=T0,T1=T1,
                     T2=T2,T01=T01[,1:N],T02=T02[,1:N],T0nu=T0nu,T10=T10[,1:N],T12=T12[,1:N],
                     T1nu=T1nu,T20=T20[,1:N],T21=T21[,1:N],T2nu=T2nu,x=x,w0=w0[,1:N],w1=w1[,1:N],
                     w2=w2[,1:N],N=N)$par
        a <- c(aop[1],aop[2],1-aop[1]-aop[2])
      }
      a <- adjust_vector(a,al,ar,zero)
      obj1_0 <- (1/N)*sum((a[1]*T0(x) + a[2]*T1(T01[,1:N]) + a[3]*T2(T02[,1:N]) - T0nu)^2 * w0[,1:N])
      obj2_0 <- (1/N)*sum((a[1]*T0(T10[,1:N]) + a[2]*T1(x) + a[3]*T2(T12[,1:N]) - T1nu)^2 * w1[,1:N])
      obj3_0 <- (1/N)*sum((a[1]*T0(T20[,1:N]) + a[2]*T1(T21[,1:N]) + a[3]*T2(x) - T2nu)^2 * w2[,1:N])
      obj_0 <- obj1_0 + obj2_0 + obj3_0
      cat(sprintf("Objective 0: Objective = %.12f\n", obj_0))
      
      result <- Main2(T0,T1,T2,a,obj_0,T10[,1:N],T12[,1:N],T1nu,w1[,1:N],T20[,1:N],T21[,1:N],T2nu,w2[,1:N],
                      T01[,1:N],T02[,1:N],T0nu,w0[,1:N],x,tol=1e-9,al,ar,zero,a_type,max_iter=500)
      
      T0_res <- result$T_res$T0_res
      T1_res <- result$T_res$T1_res
      T2_res <- result$T_res$T2_res
      a_res <- result$a_res
      RMSE_test <- sqrt((1/test)*sum(((a_res[1]*T0_res(T10[,(N+1):(N+test)])
                                       + a_res[2]*T1_res(x) + a_res[3]*T2_res(T12[,(N+1):(N+test)])
                                       - a_true[1]*T0_true(T10[,(N+1):(N+test)])-a_true[2]*T1_true(x)
                                       - a_true[3]*T2_true(T12[,(N+1):(N+test)]))^2*w1[,(N+1):(N+test)]+
                                        (a_res[1]*T0_res(T20[,(N+1):(N+test)])
                                         + a_res[2]*T1_res(T21[,(N+1):(N+test)]) + a_res[3]*T2_res(x)
                                         - a_true[1]*T0_true(T20[,(N+1):(N+test)])-a_true[2]*T1_true(T21[,(N+1):(N+test)])
                                         - a_true[3]*T2_true(x))^2*w2[,(N+1):(N+test)]+
                                        (a_res[1]*T0_res(x)
                                         + a_res[2]*T1_res(T01[,(N+1):(N+test)]) + a_res[3]*T2_res(T02[,(N+1):(N+test)])
                                         - a_true[1]*T0_true(x)-a_true[2]*T1_true(T01[,(N+1):(N+test)])
                                         - a_true[3]*T2_true(T02[,(N+1):(N+test)]))^2*w0[,(N+1):(N+test)])/3))
      
      results[[l]] <- c(sqrt(sum((a_res-a_true)^2)),sqrt(trapz(x, (T0_res(x) - T0_true(x))^2)),
                        sqrt(trapz(x, (T1_res(x) - T1_true(x))^2)),sqrt(trapz(x, (T2_res(x) - T2_true(x))^2)),
                        RMSE_test)
      
    })
    stopCluster(cl)
    
    a0_res = mean(matrix(unlist(results),nrow = 5)[1,])
    a0_sd = sd(matrix(unlist(results),nrow = 5)[1,])
    L2_T0 = mean(matrix(unlist(results),nrow = 5)[2,])
    T0_sd = sd(matrix(unlist(results),nrow = 5)[2,])
    L2_T1 = mean(matrix(unlist(results),nrow = 5)[3,])
    T1_sd = sd(matrix(unlist(results),nrow = 5)[3,])
    T2 = mean(matrix(unlist(results),nrow = 5)[4,])
    T2_sd = sd(matrix(unlist(results),nrow = 5)[4,])
    RMSE = mean(matrix(unlist(results),nrow = 5)[5,])
    RMSE_sd = sd(matrix(unlist(results),nrow = 5)[5,])
    
    
    cat(sprintf("p,a,L2_T0,L2_T1,T2:  %.3f(%.3f) %.3f(%.3f) %.3f(%.3f) %.3f(%.3f) %.3f(%.3f)\n", 
                RMSE,RMSE_sd,a0_res,a0_sd,L2_T0,T0_sd,L2_T1,T1_sd,T2,T2_sd))
    
  }
}
