## Code to reproduce figure 8
# Nonlinear confounding effects

# Import the functions needed
source("../FunctionsHDAM/FitDeconfoundedHDAM.R")
source("../FunctionsHDAM/AnalyzeFittedHDAM.R")

source("Imports.R")
library("parallel")
library(ggplot2)

f1 <- function(x){-sin(2*x)}
f2 <- function(x){2-2*tanh(x+0.5)}
f3 <- function(x){x}
f4 <- function(x){4/(exp(x)+exp(-x))}
f <- function(X){f1(X[,1])+f2(X[,2])+f3(X[,3])+f4(X[,4])}

# function to induce nonlinearity in H->X
# # linear interpolation between 0.3t and tanh(1.5t)
# nlX <- function(al, t){(1-al)*0.3*t+al*tanh(1.5*t)}

# same function as below: linear interpolation between t and abs(t)
nlX <- function(al, t){(1-al)*t+al*abs(t)}

# function to induce nonlinearity in H->Y
# linear interpolation between t and abs(t)
nlY <- function(bet, t){(1-bet)*t+bet*abs(t)}

# fix n = 400, p = 500
q <- 5
n <-  400
p <- 500

one.sim <- function(al, bet, meth, seed.val){
  set.seed(seed.val)
  # generate data
  Gamma <- matrix(runif(q*p, min = -1, max = 1), nrow=q)
  psi <- runif(q, min = 0, max = 2)
  H <- matrix(rnorm(n*q), nrow = n)
  E <- matrix(rnorm(n*p), nrow = n)
  e <- 0.5*rnorm(n)
  X <- nlX(al, H %*% Gamma) + E
  Y <- f(X) + nlY(bet, H %*% psi) + e
  lres <- FitDeconfoundedHDAM(Y, X, n.K = 5, meth = meth, cv.method = "1se", cv.k = 5, n.lambda1 = 10, n.lambda2 = 25)
  ntest <- 1000
  Htest <- matrix(rnorm(ntest * q), nrow = ntest)
  Etest <- matrix(rnorm(ntest * p), nrow = ntest)
  Xtest <- nlX(al, Htest %*% Gamma) +Etest
  fhat.res <- estimate.function(Xtest, lres)
  fXtest <- f(Xtest)
  return(c(mean((fXtest-fhat.res)^2), length(lres$active)))
}

nrep <- 50

ta <- Sys.time()
al.vec <- seq(0, 1, length.out=7)
bet.vec <- seq(0, 1, length.out=7)

l.trim <- list()
l.none <- list()

set.seed(1432)
smax <- 2100000000
seed.vec <- sample(1:smax, nrep)

for (i in 1:length(al.vec)){
  l.trim[[i]] <- list()
  l.none[[i]] <- list()
  for(l in 1:length(bet.vec)){
    al <- al.vec[i]
    bet <- bet.vec[l]
    fun = function(seed.val){return(one.sim(al, bet, "trim", seed.val))}
    l.trim[[i]][[l]] <- mclapply(seed.vec, fun, mc.cores = 60)
    fun = function(seed.val){return(one.sim(al, bet, "none", seed.val))}
    l.none[[i]][[l]] <- mclapply(seed.vec, fun, mc.cores = 60)
  }
}
te <- Sys.time()
time.needed <- te-ta

save(l.trim, l.none, time.needed, file = "ResultsNL_23_10_11.RData")

load("ResultsNL_23_10_11.RData")

trim.MSE.mean <- matrix(ncol=length(bet.vec), nrow=length(al.vec))
trim.MSE.sd <- matrix(ncol=length(bet.vec), nrow=length(al.vec))

trim.active.mean <- matrix(ncol=length(bet.vec), nrow=length(al.vec))
trim.active.sd <- matrix(ncol=length(bet.vec), nrow=length(al.vec))

none.MSE.mean <- matrix(ncol=length(bet.vec), nrow=length(al.vec))
none.MSE.sd <- matrix(ncol=length(bet.vec), nrow=length(al.vec))

none.active.mean <- matrix(ncol=length(bet.vec), nrow=length(al.vec))
none.active.sd <- matrix(ncol=length(bet.vec), nrow=length(al.vec))

for(i in 1:length(al.vec)){
  for(j in 1:length(bet.vec)){
    mses.trim <- do.call(rbind, l.trim[[i]][[j]])[, 1]
    active.trim <- do.call(rbind, l.trim[[i]][[j]])[, 2]
    mses.none <- do.call(rbind, l.none[[i]][[j]])[, 1]
    active.none <- do.call(rbind, l.none[[i]][[j]])[, 2]
    trim.MSE.mean[i,j] <- mean(mses.trim)
    trim.MSE.sd[i,j] <- sd(mses.trim)
    trim.active.mean[i,j] <- mean(active.trim)
    trim.active.sd[i,j] <- sd(active.trim)
    none.MSE.mean[i,j] <- mean(mses.none)
    none.MSE.sd[i,j] <- sd(mses.none)
    none.active.mean[i,j] <- mean(active.none)
    none.active.sd[i,j] <- sd(active.none)
  }
}


# plot MSEs
library(viridis)

par(mfrow = c(1, 2))
image(al.vec, bet.vec, trim.MSE.mean, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(0.4, 7),
      main = "Average MSE with trim transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(trim.MSE.mean[i,j] < 4){
      text(al.vec[i], bet.vec[j], round(trim.MSE.mean[i,j], 2), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(trim.MSE.mean[i,j], 2), cex = 1.5, col = "white")
    }
  
  }
}

image(al.vec, bet.vec, none.MSE.mean, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(0.4, 7),
      main = "Average MSE without transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(none.MSE.mean[i,j] < 4){
      text(al.vec[i], bet.vec[j], round(none.MSE.mean[i,j], 2), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(none.MSE.mean[i,j], 2), cex = 1.5, col = "white")
    }
    
  }
}

# plot ratio MSEtrim/MSEnone
par(mfrow=c(1,1))
ratio.MSE.mean <- trim.MSE.mean/none.MSE.mean

image(al.vec, bet.vec, ratio.MSE.mean, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),
      main = "Ratio of average MSE")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(ratio.MSE.mean[i,j] < 1){
      text(al.vec[i], bet.vec[j], round(ratio.MSE.mean[i,j], 2), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(ratio.MSE.mean[i,j], 2), cex = 1.5, col = "white")
    }
    
  }
}

# plot standard deviations
par(mfrow = c(1, 2))
image(al.vec, bet.vec, trim.MSE.sd, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(0.1, 2.6),
      main = "sd of MSEs with trim transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(trim.MSE.sd[i,j] < 1.45){
      text(al.vec[i], bet.vec[j], round(trim.MSE.sd[i,j], 2), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(trim.MSE.sd[i,j], 2), cex = 1.5, col = "white")
    }
    
  }
}

image(al.vec, bet.vec, none.MSE.sd, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(0.1, 2.6),
      main = "sd of MSEs without transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(none.MSE.sd[i,j] < 1.45){
      text(al.vec[i], bet.vec[j], round(none.MSE.sd[i,j], 2), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(none.MSE.sd[i,j], 2), cex = 1.5, col = "white")
    }
  }
}

# plot average size of active sets

par(mfrow = c(1, 2))
image(al.vec, bet.vec, trim.active.mean, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(16, 144),
      main = "Average size of active set with trim transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(trim.active.mean[i,j] < 76){
      text(al.vec[i], bet.vec[j], round(trim.active.mean[i,j], 0), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(trim.active.mean[i,j], 0), cex = 1.5, col = "white")
    }
    
  }
}

image(al.vec, bet.vec, none.active.mean, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(16, 144),
      main = "Average size of active set without transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(none.active.mean[i,j] < 76){
      text(al.vec[i], bet.vec[j], round(none.active.mean[i,j], 0), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(none.active.mean[i,j], 0), cex = 1.5, col = "white")
    }
    
  }
}
  
# plot sd of size of active sets

par(mfrow = c(1, 2))
image(al.vec, bet.vec, trim.active.sd, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(10, 37),
      main = "sd of size of active set with trim transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(trim.active.sd[i,j] < 25){
      text(al.vec[i], bet.vec[j], round(trim.active.sd[i,j], 0), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(trim.active.sd[i,j], 0), cex = 1.5, col = "white")
    }
    
  }
}

image(al.vec, bet.vec, none.active.sd, axes = FALSE, col = magma(100, direction = -1), xlab = expression(alpha), ylab = expression(beta),  zlim = c(10, 37),
      main = "sd of size of active set without transformation")

axis(1, al.vec, round(al.vec, digits = 2))
axis(2, bet.vec, round(bet.vec, digits = 2))

for (i in 1:length(al.vec)){
  for (j in 1:length(bet.vec)){
    if(none.active.sd[i,j] < 25){
      text(al.vec[i], bet.vec[j], round(none.active.sd[i,j], 0), cex = 1.5)
    }
    else{
      text(al.vec[i], bet.vec[j], round(none.active.sd[i,j], 0), cex = 1.5, col = "white")
    }
    
  }
}
