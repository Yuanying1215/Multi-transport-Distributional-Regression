# ==============================================================================
# Script: Results_figure.R
# Description: Real Data Application (Section 7).
#              Analyzes mortality data from the Human Mortality Database (HMD).
#              Compares MTDR, OT, and GOT models for forecasting male mortality 
#              distributions in 2010 based on 2005 data.
# Output: 
#   1. Prints Average Wasserstein Distance (AWD) for all models (Table 5).
#   2. Generates plots for the case study of Bulgaria (Figure 1 & 2).
#
# Data: 
#   - MortMale.RData: Male mortality data.
#   - MortFemale.RData: Female mortality data.
#
# Dependencies: 
#   - Functions.R (Ensure this file is in the same directory)
#   - Packages: pracma, fdapace, fdadensity
# ==============================================================================

library(pracma)
library(fdapace)
source("Functions.R")

load('MortMale.RData')

# ------------------------------------------------------------------------------
# Step 1: Data Pre-processing (Male)
# ------------------------------------------------------------------------------
# - Filter years: 1960-2050
# - Cut off age at 100
# - Convert densities to quantile functions (Qy: Target 2010, Qx[[1]]: Predictor 2005 Male)

n = length(year)
for (i in c(1:n)) {
  index1 <- which(year[[i]]>=1960)
  index2 <- which(year[[i]]<=2050)
  index <- intersect(index1, index2)
  year[[i]] <- year[[i]][index]
  mort[[i]] <- mort[[i]][, index]
}

## cut off at age 100
mort <- lapply(mort, function(dens){
  x <- dens[age <= 100,]
  t(t(x) / apply(x, 2, pracma::trapz,x=age[age<=100]))
})
age <- age[age <= 100]

# create qunatile functions
# qSup <- c(seq(0,0.1-0.1/50,0.1/50), seq(0.1,0.9,0.8/40), seq(0.9+0.1/50,1,0.1/50))
qSup <- seq(from=0, to = 1, length.out = 50)
qfall <- lapply(mort, function(dens) {
  apply(dens, 2, fdadensity::dens2quantile, dSup = age, qSup = qSup)
})
yearall <- year


n <- 34  # Number of countries (corresponds to sample size 'n' in paper)
p <- 2   # Number of predictors (Male 2005, Female 2005)
Qx <-  vector("list", length=(p))

x_year1 <- 2005
y_year <- 2010

Qy <- sapply(c(1:n), function(i) {
  index <- which(yearall[[i]]==y_year)
  return(qfall[[i]][,index])
})

Qx[[1]] <- sapply(c(1:n), function(i) {
  index <- which(yearall[[i]]==x_year1)
  return(qfall[[i]][,index])
})

load('MortFemale.RData')

# ------------------------------------------------------------------------------
# Step 2: Data Pre-processing (Female)
# ------------------------------------------------------------------------------
# - Apply same filtering and cutoff
# - Convert to quantile functions (Qx[[2]]: Predictor 2005 Female)

n = length(year)
for (i in c(1:n)) {
  index1 <- which(year[[i]]>=1960)
  index2 <- which(year[[i]]<=2050)
  index <- intersect(index1, index2)
  year[[i]] <- year[[i]][index]
  mort[[i]] <- mort[[i]][, index]
}

## cut off at age 100
mort <- lapply(mort, function(dens){
  x <- dens[age <= 100,]
  t(t(x) / apply(x,2,pracma::trapz,x=age[age<=100]))
})
age <- age[age <= 100]

# create quantile functions
#qSup <- c(seq(0,0.1-0.1/50,0.1/50), seq(0.1,0.9,0.8/40), seq(0.9+0.1/50,1,0.1/50))
qSup <- seq(from=0, to = 1, length.out = 50)
qfall <- lapply(mort, function(dens) {
  apply(dens, 2, fdadensity::dens2quantile, dSup = age, qSup = qSup)
})
yearall <- year

x_year2 <- 2005

Qx[[2]] <- sapply(c(1:n), function(i) {
  index <- which(yearall[[i]]==x_year2)
  return(qfall[[i]][,index])
})

# ------------------------------------------------------------------------------
# Step 3: Model Fitting & Initial Comparison (Full Dataset)
# ------------------------------------------------------------------------------
# Fits GOT, MTDR, and OT models on the entire dataset to check basic functionality.

s0=0;s1=100
m <- 101
x <- seq(s0, s1, length.out = m)

#### Method:GOT
y_mean <- rowMeans(Qy)
x_mean <-  vector("list", length=(p))
for (k in c(1:p)) {
  x_mean[[k]] <- rowMeans(Qx[[k]])
}
loss <- sapply(c(1:p), function(k){
  optimize(WassersteinLoss, c(-2, 2), Qy=Qy, y_mean=y_mean, Qx=list(Qx[[k]]), x_mean=list(x_mean[[k]]), qSup=qSup)$objective
})
Qx_old <- Qx
Qx <- Qx[order(loss)]
x_mean <- x_mean[order(loss)]
alpha <- optim(par=c(0,0), fn=WassersteinLoss, Qy=Qy, y_mean=y_mean, Qx=list(Qx[[1]], Qx[[2]]), x_mean=x_mean, qSup=qSup)$par
print(alpha)

#### Method:MTDR
bar <- y_mean
#bar <- qunif(qSup)*s1
#bar <- qbeta(qSup,0.5,0.5)*s1
#bar <- qbeta(qSup,2,3)*s1
result_mot <- MOT2(Qy,Qx_old,bar,qSup,x,tol=1e-7,inti_type=1,al=0,ar=1,zero=1e-6)
print(result_mot$a_res)

#### Method:OT
T_ot <- OT(Qy,Qx_old[[1]],qSup,x)
####

alpha_mot <- result_mot$a_res
T_res <- result_mot$T_res
k <- length(T_res)
l_mot <- rep(0,34)
l_gtr <- rep(0,34)
l_ot <- rep(0,34)
dSup <- seq(1:101)-1
for (test in 1:n) {
  X <- list(bar,Qx_old[[1]][,test], Qx_old[[2]][,test])
  mot_pred <- 0
  for (j in 1:k) {
    mot_pred <- mot_pred + alpha_mot[j]*T_res[[j]](X[[j]])
  }
  ot_pred <- T_ot(Qx_old[[1]][,test])
  gtr_pred <- WassersteinPredict(alpha, y_mean, qx=list(Qx[[1]][,test], Qx[[2]][,test]), x_mean, qSup=qSup)
  
  l_mot[test] <- sqrt(trapzRcpp(qSup, (Qy[, test]- mot_pred)^2))
  l_ot[test] <- sqrt(trapzRcpp(qSup, (Qy[, test]- ot_pred)^2))
  l_gtr[test] <- sqrt(trapzRcpp(qSup, (Qy[, test]- gtr_pred)^2))
}
mean(l_mot) # 0.4829716
mean(l_ot) # 0.598826
mean(l_gtr) # 0.5194221

# ------------------------------------------------------------------------------
# Step 4: Leave-One-Out Cross-Validation (LOOCV)
# ------------------------------------------------------------------------------
# Iterates through all 34 countries to compute the Average Wasserstein Distance (AWD)
# reported in Table 5.

### compare wasserstein distance loo
alpha_got <- matrix(0,nrow=34,ncol=2)
alpha_mot <- matrix(0,nrow = 34,ncol=3)
l_ot <- rep(0,34)
l_mot <- rep(0,34)
l_gtr <- rep(0, 34)
for (i in c(1:34)) {
  test <- c(i)
  train <- c(1:34)[-test]
  
  y_mean <- rowMeans(Qy[,train])
  x_mean <-  vector("list", length=(p))
  for (k in c(1:p)) {
    x_mean[[k]] <- rowMeans(Qx[[k]][,train])
  }
  # order predictors GOT
  loss <- sapply(c(1:p), function(k){
    optimize(WassersteinLoss, c(-2, 2), Qy=Qy[,train], y_mean=y_mean, Qx=list(Qx[[k]][,train]), x_mean=list(x_mean[[k]]), qSup=qSup)$objective
  })
  Qx_old <- Qx
  Qx <- Qx[order(loss)]
  x_mean <- x_mean[order(loss)]
  alpha <- optim(par=c(0,0), fn=WassersteinLoss, Qy=Qy[,train], y_mean=y_mean, Qx=list(Qx[[1]][,train], Qx[[2]][,train]), x_mean=x_mean, qSup=qSup)$par
  alpha_got[i,] <- alpha
  gtr_pred <- WassersteinPredict(alpha, y_mean, qx=list(Qx[[1]][,test], Qx[[2]][,test]), x_mean, qSup=qSup)
  
  #### MTDR
  bar <- y_mean
  result <- MOT2(Qy[,train],Qx=list(Qx_old[[1]][,train], Qx_old[[2]][,train]),bar,qSup,x,
                 tol=1e-7,inti_type=1,al=0,ar=1,zero=1e-6)
  #### OT
  T_ot <- OT(Qy[,train],Qx_old[[1]][,train],qSup,x)
  
  alpha_mot[i,] <- result$a_res
  T_res <- result$T_res
  k <- length(T_res)
  X <- list(bar,Qx_old[[1]][,test], Qx_old[[2]][,test])
  mot_pred <- 0
  for (j in 1:k) {
    mot_pred <- mot_pred + alpha_mot[i,j]*T_res[[j]](X[[j]])
  }
  ot_pred <- T_ot(Qx_old[[1]][,test])
  l_mot[i] <- sqrt(trapzRcpp(qSup, (Qy[, test]- mot_pred)^2))
  l_ot[i] <- sqrt(trapzRcpp(qSup, (Qy[, test]- ot_pred)^2))
  l_gtr[i] <- sqrt(trapzRcpp(qSup, (Qy[, test]- gtr_pred)^2))
}
mean(l_ot) #0.6182859
mean(l_mot) #0.5300534
mean(l_gtr) #0.5591768

# ------------------------------------------------------------------------------
# Step 5: Visualization (Figure 1 & 2)
# ------------------------------------------------------------------------------
# Case Study: Bulgaria (Index i=5)
# - Plots estimated transport maps (T0, T1, T2)
# - Plots predicted vs. actual mortality densities
############## Figure 1&2
i <- 5
country[i]
test <- c(i)
train <- c(1:34)[-test]
y_mean <- rowMeans(Qy[,train])
Qx_old <- Qx
bar <- y_mean
#bar <- qunif(qSup)*s1
result <- MOT2(Qy[,train],Qx=list(Qx_old[[1]][,train], Qx_old[[2]][,train]),bar,qSup,x,
               tol=1e-7,inti_type=1,al=0,ar=1,zero=1e-6)
print(result$a_res)
print(alpha_mot[i,])
T_res <- result$T_res
k <- length(T_res)
X <- list(bar,Qx_old[[1]][,test], Qx_old[[2]][,test])
mot_pred <- 0
for (j in 1:k) {
  mot_pred <- mot_pred + result$a_res[j]*T_res[[j]](X[[j]])
}
par(mfrow=c(1,2))
par(mai=c(0.5,0.5,0.1, 0.1))
plot(x,T_res$T1_res(x),type = "l",col="orange", lwd=2,xlab = "",ylab = "")
lines(x,x,lty=17,lwd=1)
legend("topleft", legend=c(expression(T[1]),"y=x"),
       col=c("orange","black"),bty = "n",lwd=2, lty=c(1,17), cex=1.2)
plot(0, 0, type="n",xlim=c(0, 100), ylim = c(0,0.05), xlab="", ylab="", cex.axis=1)
d1 <- frechet:::qf2pdf(Qx_old[[1]][,test], qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="brown")
d1 <- frechet:::qf2pdf(T_res$T1_res(Qx_old[[1]][,test]), qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="orange")
legend("topleft", legend=c(expression(xi[2005]^m),expression(T[1]*"#"*xi[2005]^m)),
       col=c("brown","orange"),bty = "n",lwd=2, lty=c(1,1), cex=1.2)

plot(x,T_res$T2_res(x),type = "l",col="forestgreen", lwd=2,xlab = "",ylab = "")
lines(x,x,lty=17,lwd=1)
legend("topleft", legend=c(expression(T[2]),"y=x"),
       col=c("forestgreen","black"),bty = "n",lwd=2, lty=c(1,17), cex=1.2)
plot(0, 0, type="n",xlim=c(0, 100), ylim = c(0,0.05), xlab="", ylab="", cex.axis=1)
d1 <- frechet:::qf2pdf(Qx_old[[2]][,test], qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="brown")
d1 <- frechet:::qf2pdf(T_res$T2_res(Qx_old[[2]][,test]), qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="forestgreen")
legend("topleft", legend=c(expression(xi[2005]^f),expression(T[2]*"#"*xi[2005]^f)),
       col=c("brown","forestgreen"),bty = "n",lwd=2, lty=c(1,1), cex=1.2)

plot(x,T_res$T0_res(x),type = "l",col="pink", lwd=2,xlab = "",ylab = "")
lines(x,x,lty=17,lwd=1)
legend("topleft", legend=c(expression(T[0]),"y=x"),
       col=c("pink","black"),bty = "n",lwd=2, lty=c(1,17), cex=1.2)
plot(0, 0, type="n",xlim=c(0, 100), ylim = c(0,0.05), xlab="", ylab="", cex.axis=1)
d1 <- frechet:::qf2pdf(bar, qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="brown")
d1 <- frechet:::qf2pdf(T_res$T0_res(bar), qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="pink")
legend("topleft", legend=c(expression(xi[0]),expression(T[0]*"#"*xi[0])),
       col=c("brown","pink"),bty = "n",lwd=2, lty=c(1,1), cex=1.2)

par(mfrow=c(1,1))
plot(0, 0, type="n",xlim=c(0, 100), ylim = c(0,0.05), xlab="", ylab="", cex.axis=1)
d0 <- frechet:::qf2pdf(T_res$T0_res(bar),qSup,optns = list(userBwMu = 10, outputGrid=x))
lines(d0$x, d0$y,lty=1,cex=2,col="pink", lwd=2)
d1 <- frechet:::qf2pdf(T_res$T1_res(Qx_old[[1]][,test]), qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="orange")
d1 <- frechet:::qf2pdf(T_res$T2_res(Qx_old[[2]][,test]), qSup, optns = list(userBwMu = 10, outputGrid=x))
lines(d1$x, d1$y, lty=1, lwd=2, cex=2, col="forestgreen")
d <- frechet:::qf2pdf(mot_pred, qSup, optns = list(userBwMu = 10, outputGrid=dSup))
lines(d$x, d$y, lty=1, lwd=2, cex=2, col="steelblue")
d <- frechet:::qf2pdf(Qy[, test], qSup, optns = list(userBwMu = 10, outputGrid=dSup))
lines(d$x, d$y, lty=2, lwd=2, cex=2, col="grey")
legend("topleft", legend=c(expression(widehat(eta)[2010]^m),expression(eta[2010]^m),expression(T[0]*"#"*xi[0]),expression(T[1]*"#"*xi[2005]^m)
                           ,expression(T[2]*"#"*xi[2005]^f)),
       col=c("steelblue","grey","pink", "orange","forestgreen"),bty = "n",lwd=2, lty=c(1,2,1,1), cex=1.2)









