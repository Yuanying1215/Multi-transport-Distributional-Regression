# ==============================================================================
# Script: Table4.R
# Description: Comparative simulation study for Multiple Predictors (Section 6.2).
#              Compares prediction performance (RMSE) of MTDR vs. GOT 
#              in a multi-predictor setting.
# Output: Reproduces the RMSE comparison results in Table 4.
# ==============================================================================
# --- Core Algorithm Functions ---
#' error(),mse3(): 
#' RMSE Calculations for the Proposed MTDR Model and GOT Model respectively

#' Main2():
#' Main function for MTDR

#' WassersteinLoss(),WassersteinPredict():
#' Main functions for GOT

#' parameters:
#' @param a_true
#' @param N Number of distributions/subjects. 
#'         (Note: This is denoted as 'n' in the manuscript/paper)
#' @param n Number of samples per distribution.
#'         (Note: This is denoted as 'm' in the manuscript/paper)
#' @param al,ar Lower truncation threshold/Upper truncation threshold
#' @param zero Small constant to prevent division by zero
#MTDR
error <- function(T_res,a_res,X,Y,qSup){
  #X,Y:quantile matrix:nqSup*N
  nqSup <- length(qSup)
  N <- ncol(Y)
  k <- length(T_res)
  Y_fit <- matrix(0,nrow = nqSup,ncol = N)
  mse_res <- c()
  awd_res <- c()
  for (t in 1:N) {
    for (i in 1:k) {
      Y_fit[,t] <- Y_fit[,t] + a_res[i]*T_res[[i]](X[[i]][,t])
    }
    mse_res[t] <- trapz(qSup, (Y_fit[,t]-Y[,t])^2)
    awd_res[t] <- sqrt(trapz(qSup, (Y_fit[,t]-Y[,t])^2))
    #mse[t] <- (frechet::dist4den(d1=list(x=qSup,y=Y_fit[,t]),
    #                            d2=list(x=qSup,y=Y[,t]),fctn_type = "quantile"))^2
  }
  return(list(mse = mean(mse_res), rmse = sqrt(mean(mse_res)), awd = mean(awd_res)))
}
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
lossfun_multi <- function(a,T0,T1,T2,T01,T02,T0nu,x,w0,N){
  (1/N)*sum((a[1]*T0(x) + a[2]*T1(T01) + (1-a[1]-a[2])*T2(T02) - T0nu)^2 * w0)
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
                   T2=T2_new,T01=T01,T02=T02,T0nu=T0nu,x=x,w0=w0,N=N)$par
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
#GOT
mse3 <- function(Y_fit,Y,qSup){
  N <- ncol(Y)
  mse_res <- c()
  for (t in 1:N) {
    mse_res[t] <- trapz(qSup, (Y_fit[,t]-Y[,t])^2)
  }
  return(sqrt(mean(mse_res)))
}

WassersteinLoss <- function(alpha, Qy, y_mean, Qx, x_mean, qSup){
  n <- dim(Qy)[2]
  p <- length(Qx)
  
  ## construct geodesic transports
  Tx <-  vector("list", length=(p))
  ITx <- vector("list", length=(p))
  for (k in c(1:p)) {
    Tx[[k]] <- lapply(c(1:n), function(i){splinefun(x_mean[[k]], Qx[[k]][,i], method = "natural")})
    ITx[[k]] <- lapply(c(1:n), function(i){splinefun(Qx[[k]][,i], x_mean[[k]], method = "natural")})
  }
  
  grid <- y_mean
  R = vector("list", length=(p))
  for(k in c(1:p)){
    if(alpha[p+1-k] >= 0){
      TFun <- Tx[[p+1-k]]
    }else{
      TFun <- ITx[[p+1-k]]
    }
    
    a <- floor(abs(alpha[p+1-k]))
    r <- abs(alpha[p+1-k]) - a
    
    b=1
    if(k==1){
      temp <- sapply(c(1:n), function(i){return(grid)})
    }else{
      temp <- R[[k-1]]
    }
    
    while(b<=a){
      R[[k]] <- sapply(c(1:n), function(i){return(TFun[[i]](temp[,i],0))})
      temp <- R[[k]]
      b=(b+1)
    }
    
    R[[k]] <- sapply(c(1:n), function(i){
      return(temp[,i] + r*(TFun[[i]](temp[,i],0) - temp[,i]))
    })
  }
  
  V <- R[[p]]
  loss <- mean(sapply(c(1:n), function(i){
    trapzRcpp(qSup, (V[,i]-Qy[,i])^2)
  }))#dw^2
  
  #bw <- diff(range(grid))/10
  #f <- frechet:::qf2pdf(y_mean, qSup, optns = list(userBwMu = bw, outputGrid=grid))$y
  #loss <- trapzRcpp(grid, rowMeans((V-Qy)^2)*f)
  return(loss)
}

WassersteinPredict <- function(alpha, y_mean, qx, x_mean, qSup){
  p <- length(qx)
  ## construct geodesic transports
  
  Tx <-  vector("list", length=(p))
  ITx <- vector("list", length=(p))
  
  for (k in c(1:p)) {
    Tx[[k]] <- splinefun(x_mean[[k]], qx[[k]], method = "natural")
    ITx[[k]] <- splinefun(qx[[k]], x_mean[[k]], method = "natural")
  }
  
  grid <- y_mean
  R = vector("list", length=(p))
  for(k in c(1:p)){
    if(alpha[p+1-k] >= 0){
      TFun <- Tx[[p+1-k]]
    }else{
      TFun <- ITx[[p+1-k]]
    }
    
    a <- floor(abs(alpha[p+1-k]))
    r <- abs(alpha[p+1-k]) - a
    
    if(k==1){
      temp <- grid
    }else{
      temp <- R[[k-1]]
    }
    
    b=1
    while(b<=a){
      R[[k]] <- TFun(temp,0)
      temp <- R[[k]]
      b=(b+1)
    }
    
    R[[k]] <- temp + r*(TFun(temp, 0) - temp)
  }
  
  return(R[[p]])
}


library(parallel)
library(pracma)
library(fdapace)
al <- 0.0001; ar <- 1-al
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
a_list <- list(
  c(0.3,0.35,0.35),
  c(0.2,0.2,0.6),
  c(0,0.5,0.5),
  c(0,1,0)
)
n=200
N=200
test=0.3*N
for (a_type in c(0)) {
  cat("a_type=", a_type, "\n", sep = "")
  for (a_true in a_list){
    cat("a=", a_true, "\n", sep = "")
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
      ecdf_X1 <- apply(X1q[,1:N],2,function(col){
        approxfun(col, qSup, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      eqf_X1 <- apply(X1q[,1:N],2,function(col){
        approxfun(qSup, col, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      X2 <- sapply( seq_len(N+test), function (i) {
        qbeta(runif(n), alpha2[i],beta2[i])
      })
      X2q <- sapply(seq_len(N+test), function (i) {
        quantile(X2[,i],qSup)
      })
      ecdf_X2 <- apply(X2q[,1:N],2,function(col){
        approxfun(col, qSup, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      eqf_X2 <- apply(X2q[,1:N],2,function(col){
        approxfun(qSup, col, method = "linear",
                  yleft = 0, yright = 1, ties = "ordered")
      })
      
      ecdf_bar <- function(x){x}
      eqf_bar <- function(x){x}
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
      
      T10 <- sapply(seq_len(N), function(i){
        eqf_bar(ecdf_X1[[i]](x))
      })
      T12 <- sapply(seq_len(N), function(i){
        eqf_X2[[i]](ecdf_X1[[i]](x))
      })
      T1nu <- sapply(seq_len(N), function(i){
        eqf_Y[[i]](ecdf_X1[[i]](x))
      })
      w1 <- sapply(seq_len(N), function(i){
        ecdf_X1[[i]](x)
      })
      w1 <- w1 - rbind(0,w1[-nrow(w1),])
      T20 <- sapply(seq_len(N), function(i){
        eqf_bar(ecdf_X2[[i]](x))
      })
      T21 <- sapply(seq_len(N), function(i){
        eqf_X1[[i]](ecdf_X2[[i]](x))
      })
      T2nu <- sapply(seq_len(N), function(i){
        eqf_Y[[i]](ecdf_X2[[i]](x))
      })
      w2 <- sapply(seq_len(N), function(i){
        ecdf_X2[[i]](x)
      })
      w2 <- w2 - rbind(0,w2[-nrow(w2),])
      T01 <- sapply(seq_len(N), function(i){
        eqf_X1[[i]](ecdf_bar(x))
      })
      T02 <- sapply(seq_len(N), function(i){
        eqf_X2[[i]](ecdf_bar(x))
      })
      T0nu <- sapply(seq_len(N), function(i){
        eqf_Y[[i]](ecdf_bar(x))
      })
      w0 <- sapply(seq_len(N), function(i){
        ecdf_bar(x)
      })
      w0 <- w0 - rbind(0,w0[-nrow(w0),])
      
      T2 <- approxfun(x,x,method="linear")
      T1 <- approxfun(x,x,method="linear")
      T0 <- approxfun(x,x,method="linear")
      #####
      
      if (a_type==0){
        a <- alpha_l(list(T0,T1,T2),list(T01,T02,T0nu),x,w0)
      }else{
        aop <- optim(par=c(0,0),fn=lossfun_multi,T0=T0,T1=T1,
                     T2=T2,T01=T01,T02=T02,T0nu=T0nu,x=x,w0=w0,N=N)$par
        a <- c(aop[1],aop[2],1-aop[1]-aop[2])
      }
      a <- adjust_vector(a,al,ar,zero)
      obj1_0 <- (1/N)*sum((a[1]*T0(x) + a[2]*T1(T01) + a[3]*T2(T02) - T0nu)^2 * w0)
      obj2_0 <- (1/N)*sum((a[1]*T0(T10) + a[2]*T1(x) + a[3]*T2(T12) - T1nu)^2 * w1)
      obj3_0 <- (1/N)*sum((a[1]*T0(T20) + a[2]*T1(T21) + a[3]*T2(x) - T2nu)^2 * w2)
      obj_0 <- obj1_0 + obj2_0 + obj3_0
      cat(sprintf("Objective 0: Objective = %.12f\n", obj_0))
      
      result <- Main2(T0,T1,T2,a,obj_0,T10,T12,T1nu,w1,T20,T21,T2nu,w2,
                      T01,T02,T0nu,w0,x,tol=1e-9,al,ar,zero,a_type,max_iter=500)
      
      T_res <- result$T_res
      a_res <- result$a_res
      
      #GOT 
      y_mean <- rowMeans(Y[,1:N])
      x_mean <-  vector("list", length=(2))
      Qx <- list(X1q,X2q)
      for (k in c(1:2)) {
        x_mean[[k]] <- rowMeans(Qx[[k]][,1:N])
      }
      loss <- sapply(c(1:2), function(k){
        optimize(WassersteinLoss, c(-2, 2), Qy=Y[,1:N], y_mean=y_mean, Qx=list(Qx[[k]][,1:N]), x_mean=list(x_mean[[k]]), qSup=qSup)$objective
      })
      Qx_old <- Qx
      Qx <- Qx[order(loss)]
      #Qx_old <- Qx_old[rev(order(loss))]
      x_mean_old <- x_mean
      x_mean <- x_mean[order(loss)]
      #x_mean_old <- x_mean_old[rev(order(loss))]
      alpha_got <- optim(par=c(0,0), fn=WassersteinLoss, Qy=Y[,1:N], y_mean=y_mean, Qx=list(Qx[[1]][,1:N], Qx[[2]][,1:N]), x_mean=x_mean, qSup=qSup)
      
      gtr_pred <- sapply(1:test,function(i)
      {WassersteinPredict(alpha_got$par, y_mean=y_mean, qx=list(Qx[[1]][,i+N], Qx[[2]][,i+N]),
                          x_mean=x_mean, qSup=qSup)})
      
      loss1 <- result$loss_res
      loss2 <- alpha_got$value
      a_GOT <- alpha_got$par
      
      MSE1 <- error(T_res,a_res,
                    list(matrix(rep(qSup,test),ncol=test),
                         X1q[,(N+1):(N+test)],X2q[,(N+1):(N+test)]),Y[,(N+1):(N+test)],qSup)
      MSE2 <- mse3(gtr_pred,Y[,(N+1):(N+test)],qSup)
      
      results[[l]] <- c(MSE1$rmse,MSE2)
    })
    stopCluster(cl)
    
    MSE1 = mean(matrix(unlist(results),nrow = 2)[1,])
    MSE1_sd = sd(matrix(unlist(results),nrow = 2)[1,])
    MSE2 = mean(matrix(unlist(results),nrow = 2)[2,])
    MSE2_sd = sd(matrix(unlist(results),nrow = 2)[2,])
    
    cat(sprintf("RMSE1,RMSE2: %.3f(%.3f)  %.3f(%.3f) \n",
                MSE1,MSE1_sd,MSE2,MSE2_sd))
  }
}


