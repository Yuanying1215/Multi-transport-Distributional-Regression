# ==============================================================================
# Script: TableS1.R
# Description: Supplementary simulation study (Single Predictor).
#              Compares the proposed MTDR model against the GOT (Geodesic 
#              Optimal Transport) model and OT model under the data-generating process 
#              specifc to the GOT framework.
# Output: Reproduces the RMSE comparison results reported in Supplementary 
#         Material Table S1.
#
# Scenario: Data generated according to the GOT model structure.
# ==============================================================================

# --- Core Algorithm Functions ---
#' mse(),mse2(),mse3(): 
#' RMSE Calculations for the Proposed MTDR Model, OT Model and GOT Model respectively

#' MOT():
#' Main function for MTDR

#' OT(): 
#' Main function for OT

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
mse <- function(T_res,a_res,X,Y,qSup){
  #X,Y:quantile matrix:nqSup*N
  nqSup <- length(qSup)
  N <- ncol(Y)
  k <- length(T_res)
  Y_fit <- matrix(0,nrow = nqSup,ncol = N)
  mse_res <- c()
  for (t in 1:N) {
    for (i in 1:k) {
      Y_fit[,t] <- Y_fit[,t] + a_res[i]*T_res[[i]](X[[i]][,t])
    }
    mse_res[t] <- trapz(qSup, (Y_fit[,t]-Y[,t])^2)
    #mse[t] <- (frechet::dist4den(d1=list(x=qSup,y=Y_fit[,t]),
    #                            d2=list(x=qSup,y=Y[,t]),fctn_type = "quantile"))^2
  }
  return(sqrt(mean(mse_res)))
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
lossfun <- function(a0,T0,T1,T01,T0nu,x,w0,N){
  (1/N)*sum((a0*T0(x) + (1-a0)*T1(T01) - T0nu)^2 * w0)
}
Main_new <- function(T0,T1,a0,obj_0,T10,T1nu,w1,T01,T0nu,w0,x,tol,
                     al,ar,zero,a_type,max_iter=5000)
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
    } else {a <- optimize(lossfun,c(0,1),T0_new,T1_new,T01,
                          T0nu,x,w0,N)$minimum
    }
    if (a < al){a <- zero}
    if (a > ar){a <- 1 - zero}
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
MOT <- function(Qy,Qx,bar,qSup,x,tol,inti_type,al,ar,zero,a_type){
  m <- length(x)
  s0 <- min(x);s1 <- max(x)
  N <- ncol(Qy)
  ecdf_bar <- approxfun(bar, qSup, method = "linear",
                        yleft = 0, yright = 1, ties = "ordered")
  eqf_bar <-  approxfun(qSup, bar, method = "linear",
                        yleft = s0, yright = s1, ties = "ordered")
  #   ecdf_bar <- function(x) {x}
  #   eqf_bar <- function(x) {x}
  ecdf_Y <- apply(Qy,2,function(col){
    approxfun(col, qSup, method = "linear", 
              yleft = 0, yright = 1, ties = "ordered")
  })
  eqf_Y <- apply(Qy,2,function(col){
    approxfun(qSup, col, method = "linear", 
              yleft = s0, yright = s1, ties = "ordered")
  })
  ecdf_X1 <- apply(Qx,2,function(col){
    approxfun(col, qSup, method = "linear", 
              yleft = 0, yright = 1, ties = "ordered")
  })
  eqf_X1 <- apply(Qx,2,function(col){
    approxfun(qSup, col, method = "linear", 
              yleft = s0, yright = s1, ties = "ordered")
  })
  T10 <- sapply(seq_len(N), function(i){
    eqf_bar(ecdf_X1[[i]](x))
  })
  T1nu <- sapply(seq_len(N), function(i){
    eqf_Y[[i]](ecdf_X1[[i]](x))
  })
  w1 <- sapply(seq_len(N), function(i){
    ecdf_X1[[i]](x)
  })
  w1 <- w1 - rbind(0,w1[-nrow(w1),])
  T01 <- sapply(seq_len(N), function(i){
    eqf_X1[[i]](ecdf_bar(x))
  })
  T0nu <- sapply(seq_len(N), function(i){
    eqf_Y[[i]](ecdf_bar(x))
  })
  w0 <- sapply(seq_len(N), function(i){
    (ecdf_bar(x))
  })
  w0 <- w0 - rbind(0,w0[-nrow(w0),])
  if (inti_type==0){
    T1 <- approxfun(x, c(s0,sort(runif((m-2),s0,s1)),s1), method = "linear")
    T0 <- approxfun(x, c(s0,sort(runif((m-2),s0,s1)),s1), method = "linear")
  }else{
    T0 <- approxfun(x, x, method = "linear")
    T1 <- approxfun(x, x, method = "linear")
  }
  if (a_type == 0){
    a0 <- sum(rowMeans((T0nu - T1(T01))*(T0(x)-T1(T01)))*w0[,1])/
      sum(rowMeans((matrix((T0(x)- T1(T01)),nrow=m))^2)*w0[,1])
  } else {a0 <- optimize(lossfun,c(0,1),T0,T1,T01,
                         T0nu,x,w0,N)$minimum
  }
  if (a0 < al){a0 <- zero}
  if (a0 > ar){a0 <- 1 - zero}
  obj1_0 <- (1/N)*sum((a0*T0(x) + (1-a0)*T1(T01) - T0nu)^2 * w0)
  obj2_0 <- (1/N)*sum((a0*T0(T10) + (1-a0)*T1(x) - T1nu)^2 * w1)
  obj_0 <- obj1_0 + obj2_0
  cat(sprintf("Objective 0: Objective = %.12f\n", obj_0))
  result <- Main_new(T0,T1,a0,obj_0,T10,T1nu,w1,T01,T0nu,w0,x,tol,
                     al,ar,zero,a_type)
  return(result)
}

#OT
mse2 <- function(T_res,X,Y,qSup){
  #X,Y:quantile matrix:nqSup*N
  nqSup <- length(qSup)
  N <- ncol(Y)
  Y_fit <- matrix(0,nrow = nqSup,ncol = N)
  mse_res <- c()
  for (t in 1:N) {
    Y_fit[,t] <- T_res(X[,t])
    mse_res[t] <- trapz(qSup, (Y_fit[,t]-Y[,t])^2)
  }
  return(sqrt(mean(mse_res)))
}
OT <- function(Qy,qx,qSup,x){
  s0 <- min(x);s1 <- max(x)
  N <- ncol(Qy)
  ecdf_Y <- apply(Qy,2,function(col){
    approxfun(col, qSup, method = "linear", 
              yleft = 0, yright = 1, ties = "ordered")
  })
  eqf_Y <- apply(Qy,2,function(col){
    approxfun(qSup, col, method = "linear", 
              yleft = s0, yright = s1, ties = "ordered")
  })
  ecdf_X1 <- apply(qx,2,function(col){
    approxfun(col, qSup, method = "linear", 
              yleft = 0, yright = 1, ties = "ordered")
  })
  eqf_X1 <- apply(qx,2,function(col){
    approxfun(qSup, col, method = "linear", 
              yleft = s0, yright = s1, ties = "ordered")
  })
  T1nu <- sapply(seq_len(N), function(i){
    eqf_Y[[i]](ecdf_X1[[i]](x))
  })
  w1 <- sapply(seq_len(N), function(i){
    ecdf_X1[[i]](x)
  })
  w1 <- w1 - rbind(0,w1[-nrow(w1),])
  T2_grid <- gpava_new(x,T1nu,weights = w1,s0,s1)$x 
  T2_res <- approxfun(x, T2_grid, method = "linear")
  return(T2_res)
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
library(pracma)#
library(fdapace)#
n = 200
N = 200
a_type = 1;inti_type=1
al=0.0001;ar <- 1-al
zero <- 1e-6
zeta <- function(x, k) {
  if (k == 0) {
    return(x)
  }
  x - sin(pi * k * x) / (pi * abs(k))
}

nqSup <- 101
qSup <- seq(0, 1, length.out = nqSup)
m <- 101
x <- seq(0, 1, length.out = m)
test <- 0.3*N
B <-  100
seeds <- seq(1,B,1)
for (alpha in c(0,0.5,1)){
  print(alpha)
  results <- list()
  cl <- makeCluster(20,type="FORK")  
  results <- parLapply(cl,1:B,function(l){
    set.seed(seeds[l])
    X_mean <- sapply(seq_len(N+test), function(i){
      qbeta(runif(n), 2,4)#
    })
    #X_mean <- qbeta(runif(n), 2,4)
    #K <- sample(-6:6,N+test,replace = T)
    K <- sample(c(-3,-2,0,2,3),N+test,replace=T)
    X_sample <- sapply(seq_len(N+test), function(i){
      zeta(X_mean[,i],K[i])
    })
    
    Y_mean <- sapply(seq_len(N+test), function(i){
      qbeta(runif(n), 6,2)#
    })
    if (alpha >= 0){
      Y_res <- sapply(seq_len(N+test), function(i){
        (1-alpha)*Y_mean[,i] + alpha*zeta(Y_mean[,i],K[i])
      })
    }else{
      Y_res <- sapply(seq_len(N+test), function(i){
        (1+alpha)*Y_mean[,i] - alpha*qbeta(zeta(qbeta(Y_mean[,i],2,4),K[i]),2,4)
      })
    }
    
    ###
    #E <- sample(-50:50,N+test,replace = T)
    values <- setdiff(-20:20, c(0))
    E <- sample(values,N+test,replace = T)
    Y_sample <- sapply(seq_len(N+test), function(i){
      zeta(Y_res[,i],E[i])
    })
    
    #GOT
    X1q <- sapply(seq_len(N+test), function (i) {
      quantile(X_sample[,i],qSup)
    })
    x_mean <- list(rowMeans(X1q[,1:N]))
    Yq <- sapply(seq_len(N+test), function (i) {
      quantile(Y_sample[,i],qSup)
    })
    y_mean <- rowMeans(Yq[,1:N])
    alpha_res <- optimize(WassersteinLoss, c(-2, 2), Qy=Yq[,1:N], y_mean=y_mean, Qx=list(X1q[,1:N]), x_mean=x_mean, qSup=qSup)
    alp <- alpha_res$minimum
    #MTDR
    result <- MOT(Yq[,1:N],X1q[,1:N],qSup,qSup,x,tol=1e-9,inti_type,al,ar,zero,a_type)
    T0_res <- result$T0_res
    T1_res <- result$T1_res
    a_res <- result$a_res
    #OT
    ecdf_Y <- apply(Yq[,1:N],2,function(col){
      approxfun(col, qSup, method = "linear", 
                yleft = 0, yright = 1, ties = "ordered")
    })
    eqf_Y <- apply(Yq[,1:N],2,function(col){
      approxfun(qSup, col, method = "linear", 
                yleft = 0, yright = 1, ties = "ordered")
    })
    ecdf_X1 <- apply(X1q[,1:N],2,function(col){
      approxfun(col, qSup, method = "linear", 
                yleft = 0, yright = 1, ties = "ordered")
    })
    eqf_X1 <- apply(X1q[,1:N],2,function(col){
      approxfun(qSup, col, method = "linear", 
                yleft = 0, yright = 1, ties = "ordered")
    })
    T1nu <- sapply(seq_len(N), function(i){
      eqf_Y[[i]](ecdf_X1[[i]](x))
    })
    w1 <- sapply(seq_len(N), function(i){
      ecdf_X1[[i]](x)
    })
    w1 <- w1 - rbind(0,w1[-nrow(w1),])
    T2_grid <- gpava_new(x,T1nu,weights = w1,0,1)$x 
    T2_res <- approxfun(x, T2_grid, method = "linear")
    
    gtr_pred <- sapply(1:test,function(i)
    {WassersteinPredict(alp, y_mean=y_mean, qx=list(X1q[,i+N]),
                        x_mean=x_mean, qSup=qSup)})
    a_MOT <- a_res
    a_GOT <- alp
    loss1 <- result$loss_res
    loss2 <- (1/N)*sum((T2_res(x) - T1nu)^2 * w1)
    loss3 <- alpha_res$objective
    MSE1 <- mse(list(T0_res,T1_res),c(a_res,1-a_res),
                list(matrix(rep(qSup,test),ncol=test),
                     X1q[,(N+1):(N+test)]),Yq[,(N+1):(N+test)],qSup)
    MSE2 <- mse2(T2_res,X1q[,(N+1):(N+test)],Yq[,(N+1):(N+test)],qSup)
    
    MSE3 <- mse3(gtr_pred,Yq[,(N+1):(N+test)],qSup)
    results[[l]] <- c(loss1,loss2,loss3,MSE1,MSE2,MSE3,a_MOT,a_GOT)
  })
  stopCluster(cl)
  loss1 = mean(matrix(unlist(results),nrow = 8)[1,])
  loss1_sd = sd(matrix(unlist(results),nrow = 8)[1,])
  loss2 = mean(matrix(unlist(results),nrow = 8)[2,])
  loss2_sd = sd(matrix(unlist(results),nrow = 8)[2,])
  loss3 = mean(matrix(unlist(results),nrow = 8)[3,])
  loss3_sd = sd(matrix(unlist(results),nrow = 8)[3,])
  MSE1 = mean(matrix(unlist(results),nrow = 8)[4,])
  MSE1_sd = sd(matrix(unlist(results),nrow = 8)[4,])
  MSE2 = mean(matrix(unlist(results),nrow = 8)[5,])
  MSE2_sd = sd(matrix(unlist(results),nrow = 8)[5,])
  MSE3 = mean(matrix(unlist(results),nrow = 8)[6,])
  MSE3_sd = sd(matrix(unlist(results),nrow = 8)[6,])
  aMOT = mean(matrix(unlist(results),nrow = 8)[7,])
  aMOT_sd = sd(matrix(unlist(results),nrow = 8)[7,])
  aGOT = mean(matrix(unlist(results),nrow = 8)[8,])
  aGOT_sd = sd(matrix(unlist(results),nrow = 8)[8,])
  # cat(sprintf("loss1,loss2,loss3,MSE1,MSE2,MSE3: %.4f(%.4f) %.4f(%.4f) 
  #             %.4f(%.4f) %.4f(%.4f) %.4f(%.4f) %.4f(%.4f)\n", 
  #            loss1,loss1_sd,loss2,loss2_sd,loss3,loss3_sd,
  #             MSE1,MSE1_sd,MSE2,MSE2_sd,MSE3,MSE3_sd))
  # cat(sprintf("alpha_MOT,alpha_GOT: %.4f(%.4f) %.4f(%.4f)\n", 
  #             aMOT,aMOT_sd,
  #             aGOT,aGOT_sd))
  cat(sprintf("MSE1,MSE2,MSE3: %.3f(%.3f) %.3f(%.3f) %.3f(%.3f)\n", 
              MSE1,MSE1_sd,MSE2,MSE2_sd,MSE3,MSE3_sd))
}
