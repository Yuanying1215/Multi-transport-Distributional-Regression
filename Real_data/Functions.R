# ==============================================================================
# Script: Functions.R
# Description: Helper functions for the Real Data Application (Section 7).
#              Contains core algorithms for MTDR (gPAVA, Algorithm 1) and 
#              comparative methods (OT, GOT).
# Usage: This script is sourced by 'Results_figure.R'.
# ==============================================================================

# --- Core Algorithm Functions ---
#' gpava_new():
#' Generalized Pool Adjacent Violators Algorithm (gPAVA)
#' Solves the weighted isotonic regression sub-problem for updating transport maps.

#' lossfun_multi(),lossfun():
#' Loss Functions
#' Calculates the objective function value for optimization.

#' Main2(),MOT2(),Main_new(),MOT():
#' Main Iterative Optimization (Algorithm 1)
#' Alternates between updating transport maps (via gPAVA) and weights.

#' OT():
#' Function for OT Model

#' WassersteinLoss(),WassersteinPredict():
#' Main functions for GOT Model

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
lossfun_multi <- function(a,T0,T1,T2,T01,T02,T0nu,T10,T12,T1nu,T20,T21,T2nu,x,w0,w1,w2,N){
  (1/N)*sum((a[1]*T0(x) + a[2]*T1(T01) + (1-a[1]-a[2])*T2(T02) - T0nu)^2 * w0)+
    (1/N)*sum((a[1]*T0(T10) + a[2]*T1(x) + (1-a[1]-a[2])*T2(T12) - T1nu)^2 * w1)+
    (1/N)*sum((a[1]*T0(T20) + a[2]*T1(T21) + (1-a[1]-a[2])*T2(x) - T2nu)^2 * w2)
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
Main2 <- function(T0,T1,T2,a,obj_0,T10,T12,T1nu,w1,T20,T21,T2nu,w2,
                  T01,T02,T0nu,w0,x,tol,al,ar,zero,max_iter=500)
{
  a_type = 1
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
      aop <- optim(par=c(0.35,0.35),fn=lossfun_multi,T0=T0_new,T1=T1_new,
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
MOT2 <- function(Qy,Qx,bar,qSup,x,tol,inti_type,al,ar,zero){
  a_type <- 1
  m <- length(x)
  s0 <- min(x);s1 <- max(x)
  N <- ncol(Qy)
  ecdf_bar <- approxfun(bar, qSup, method = "linear",
                        yleft = 0, yright = 1, ties = "ordered")
  eqf_bar <-  approxfun(qSup, bar, method = "linear",
                        yleft = s0, yright = s1, ties = "ordered")
  ecdf_Y <- apply(Qy,2,function(col){
    approxfun(col, qSup, method = "linear", 
              yleft = 0, yright = 1, ties = "ordered")
  })
  eqf_Y <- apply(Qy,2,function(col){
    approxfun(qSup, col, method = "linear", 
              yleft = s0, yright = s1, ties = "ordered")
  })
  ecdf_X1 <- apply(Qx[[1]],2,function(col){
    approxfun(col, qSup, method = "linear", 
              yleft = 0, yright = 1, ties = "ordered")
  })
  eqf_X1 <- apply(Qx[[1]],2,function(col){
    approxfun(qSup, col, method = "linear", 
              yleft = s0, yright = s1, ties = "ordered")
  })
  ecdf_X2 <- apply(Qx[[2]],2,function(col){
    approxfun(col, qSup, method = "linear", 
              yleft = 0, yright = 1, ties = "ordered")
  })
  eqf_X2 <- apply(Qx[[2]],2,function(col){
    approxfun(qSup, col, method = "linear", 
              yleft = s0, yright = s1, ties = "ordered")
  })
  T12 <- sapply(seq_len(N), function(i){
    eqf_X2[[i]](ecdf_X1[[i]](x))
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
  T02 <- sapply(seq_len(N), function(i){
    eqf_X2[[i]](ecdf_bar(x))
  })
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
  T21 <- sapply(seq_len(N), function(i){
    eqf_X1[[i]](ecdf_X2[[i]](x))
  })
  T20 <- sapply(seq_len(N), function(i){
    eqf_bar(ecdf_X2[[i]](x))
  })
  T2nu <- sapply(seq_len(N), function(i){
    eqf_Y[[i]](ecdf_X2[[i]](x))
  })
  w2 <- sapply(seq_len(N), function(i){
    (ecdf_X2[[i]](x))
  })
  w2 <- w2 - rbind(0,w2[-nrow(w2),])
  if (inti_type==0){
    T1 <- approxfun(x, c(s0,sort(runif((m-2),s0,s1)),s1), method = "linear")
    T0 <- approxfun(x, c(s0,sort(runif((m-2),s0,s1)),s1), method = "linear")
    T2 <- approxfun(x, c(s0,sort(runif((m-2),s0,s1)),s1), method = "linear")
  }else{
    T0 <- approxfun(x, x, method = "linear")
    T1 <- approxfun(x, x, method = "linear")
    T2 <- approxfun(x, x, method = "linear")
  }
  if (a_type==0){
    a <- alpha_l(list(T0,T1,T2),list(T01,T02,T0nu),x,w0)
  }else{
    ao <- optim(par=c(0.35,0.35),fn=lossfun_multi,T0=T0,T1=T1,
                T2=T2,T01=T01,T02=T02,T0nu=T0nu,T10=T10,T12=T12,
                T1nu=T1nu,T20=T20,T21=T21,T2nu=T2nu,x=x,w0=w0,w1=w1,
                w2=w2,N=N)$par
    a <- matrix(c(ao[1],ao[2],1-ao[1]-ao[2]),nrow=3)
  }
  a <- adjust_vector(a,al,ar,zero)
  obj1_0 <- (1/N)*sum((a[1,1]*T0(x) + a[2,1]*T1(T01) + a[3,1]*T2(T02) - T0nu)^2 * w0)
  obj2_0 <- (1/N)*sum((a[1,1]*T0(T10) + a[2,1]*T1(x) + a[3,1]*T2(T12) - T1nu)^2 * w1)
  obj3_0 <- (1/N)*sum((a[1,1]*T0(T20) + a[2,1]*T1(T21) + a[3,1]*T2(x) - T2nu)^2 * w2)
  obj_0 <- obj1_0 + obj2_0 + obj3_0
  cat(sprintf("Objective 0: Objective = %.12f\n", obj_0))
  result <- Main2(T0,T1,T2,a,obj_0,T10,T12,T1nu,w1,T20,T21,T2nu,w2,
                  T01,T02,T0nu,w0,x,tol,al,ar,zero)
  return(result)
}

lossfun <- function(a0,T0,T1,T01,T0nu,w0,T10,T1nu,x,w1,N){
  (1/N)*sum((a0*T0(T10) + (1-a0)*T1(x) - T1nu)^2 * w1)+
    (1/N)*sum((a0*T0(x) + (1-a0)*T1(T01) - T0nu)^2 * w0)
}
Main_new <- function(T0,T1,a0,obj_0,T10,T1nu,w1,T01,T0nu,w0,x,tol,
                     al,ar,zero,max_iter=500)
{
  a_type <- 1
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
    } else {a <- optimize(lossfun,c(0,1),T0_new,T1_new,
                          T01,T0nu,w0,T10,T1nu,x,w1,N)$minimum
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
MOT <- function(Qy,Qx,bar,qSup,x,tol,inti_type,al,ar,zero){
  a_type <- 1
  m <- length(x)
  s0 <- min(x);s1 <- max(x)
  N <- ncol(Qy)
  ecdf_bar <- approxfun(bar, qSup, method = "linear",
                        yleft = 0, yright = 1, ties = "ordered")
  eqf_bar <-  approxfun(qSup, bar, method = "linear",
                        yleft = s0, yright = s1, ties = "ordered")
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
    a0 <- sum((T0nu - T1(T01))*(T0(x)-T1(T01))*w0)/
      sum(rowMeans((matrix((T0(x)- T1(T01)),nrow=m))^2)*w0)
    if (a0 < al){a0 <- zero}
    if (a0 > ar){a0 <- 1 - zero}
  } else {a0 <- optimize(lossfun,c(0,1),T0,T1,T01,T0nu,w0,T10,T1nu,x,w1,N)$minimum
  }
  obj1_0 <- (1/N)*sum((a0*T0(x) + (1-a0)*T1(T01) - T0nu)^2 * w0)
  obj2_0 <- (1/N)*sum((a0*T0(T10) + (1-a0)*T1(x) - T1nu)^2 * w1)
  obj_0 <- obj1_0 + obj2_0
  cat(sprintf("Objective 0: Objective = %.12f\n", obj_0))
  result <- Main_new(T0,T1,a0,obj_0,T10,T1nu,w1,T01,T0nu,w0,x,tol,
                     al,ar,zero)
  return(result)
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
