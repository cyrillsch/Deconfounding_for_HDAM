## Code to reproduce Figures 2 and 3
# Vary the sample size n

# Import the functions needed
source("../FunctionsHDAM/FitDeconfoundedHDAM.R")
source("../FunctionsHDAM/AnalyzeFittedHDAM.R")

source("Imports.R")
library("parallel")
library(ggplot2)
detectCores()

f1 <- function(x){-sin(2*x)}
f2 <- function(x){2-2*tanh(x+0.5)}
f3 <- function(x){x}
f4 <- function(x){4/(exp(x)+exp(-x))}
f <- function(X){f1(X[,1])+f2(X[,2])+f3(X[,3])+f4(X[,4])}

# fix p = 300, q = 5
q <- 5
p <- 300

one.sim <- function(n, meth, rho = NULL, seed.val){
  # generate data
  set.seed(seed.val)
  Gamma <- matrix(runif(q*p, min = -1, max = 1), nrow=q)
  psi <- runif(q, min = 0, max = 2)
  H <- matrix(rnorm(n*q), nrow = n)
  # E should be Toeplitz if rho is specified
  if(is.null(rho)){
    E <- matrix(rnorm(n*p), nrow = n)
  } else {
    E <- matrix(rnorm(n*p), nrow = n)
    Toe <- toeplitz(rho^(0:(p-1)))
    # t(R.toe)%*%R.toe = Toe
    R.toe <- chol(Toe)
    # E has rows i.i.d with covariance Toe
    E <- E %*% R.toe
  }
  e <- 0.5*rnorm(n)
  X <- H %*% Gamma + E
  Y <- f(X) + H %*% psi + e
  lres <- FitDeconfoundedHDAM(Y, X, n.K = 5, meth = meth, cv.method = "1se", cv.k = 5, n.lambda1 = 10, n.lambda2 = 25)
  ntest <- 1000
  Htest <- matrix(rnorm(ntest * q), nrow = ntest)
  Etest <- matrix(rnorm(ntest * p), nrow = ntest)
  Xtest <- Htest %*% Gamma +Etest
  fhat.res <- estimate.function(Xtest, lres)
  fXtest <- f(Xtest)
  return(c(mean((fXtest-fhat.res)^2), length(lres$active)))
}

nrep <- 100

ta <- Sys.time()
n.vec <- c(50, 100, 200, 400, 800)

l.VarN.trim.NULL <- list()
l.VarN.none.NULL <- list()
l.VarN.trim.rho04 <- list()
l.VarN.none.rho04 <- list()
l.VarN.trim.rho08 <- list()
l.VarN.none.rho08 <- list()

set.seed(1432)
smax <- 2100000000
seed.vec <- sample(1:smax, nrep)

for (i in 1:length(n.vec)){
  n <- n.vec[[i]]
  rho = NULL
  fun = function(seed.val){return(one.sim(n, "trim", rho = NULL, seed.val = seed.val))}
  l.VarN.trim.NULL[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  fun = function(seed.val){return(one.sim(n, "none", rho = NULL, seed.val = seed.val))}
  l.VarN.none.NULL[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  # rho = 0.4
  fun = function(seed.val){return(one.sim(n, "trim", rho = 0.4, seed.val = seed.val))}
  l.VarN.trim.rho04[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  fun = function(seed.val){return(one.sim(n, "none", rho = 0.4, seed.val = seed.val))}
  l.VarN.none.rho04[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  # rho = 0.8
  fun = function(seed.val){return(one.sim(n, "trim", rho = 0.8, seed.val = seed.val))}
  l.VarN.trim.rho08[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  fun = function(seed.val){return(one.sim(n, "none", rho = 0.8, seed.val = seed.val))}
  l.VarN.none.rho08[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
}
te <- Sys.time()
time.needed <- te-ta

save(l.VarN.trim.NULL, l.VarN.none.NULL, l.VarN.trim.rho04, l.VarN.none.rho04, l.VarN.trim.rho08, l.VarN.none.rho08, time.needed, file = "VarN_2023_10_12.RData")

# Plot results
load("VarN_2023_10_12.RData")

# Components of E independent
resTrim.NULL <- data.frame(matrix(ncol = 4, nrow = nrep * length(n.vec)))
colnames(resTrim.NULL) <- c("n", "MSE", "s.active", "meth")

for(j in 1: length(n.vec)){
  for(i in 1:nrep){
    resTrim.NULL$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim.NULL$n[(j-1)*nrep + i] <- n.vec[j]
    resTrim.NULL$MSE[(j-1)*nrep + i] <- l.VarN.trim.NULL[[j]][[i]][1]
    resTrim.NULL$s.active[(j-1)*nrep + i] <- l.VarN.trim.NULL[[j]][[i]][2]
  }
}

resNone.NULL <- data.frame(matrix(ncol = 4, nrow = nrep * length(n.vec)))
colnames(resNone.NULL) <- c("n", "MSE", "s.active", "meth")

for(j in 1: length(n.vec)){
  for(i in 1:nrep){
    resNone.NULL$meth[(j-1)*nrep + i] <- "naive"
    resNone.NULL$n[(j-1)*nrep + i] <- n.vec[j]
    resNone.NULL$MSE[(j-1)*nrep + i] <- l.VarN.none.NULL[[j]][[i]][1]
    resNone.NULL$s.active[(j-1)*nrep + i] <- l.VarN.none.NULL[[j]][[i]][2]
  }
}

resTrim.NULL$n <- factor(resTrim.NULL$n, levels = n.vec)
resTrim.NULL$meth <- as.factor(resTrim.NULL$meth)
resNone.NULL$n <- factor(resNone.NULL$n, levels = n.vec)
resNone.NULL$meth <- as.factor(resNone.NULL$meth)

resTot.NULL <- rbind(resTrim.NULL, resNone.NULL)

library(gridExtra)
p <- ggplot(resTot.NULL, aes(x=n, y=MSE, fill=meth))+geom_violin(scale = "width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("n") + ylab("MSE") + ggtitle("MSE of f with p=300, q=5, s=4, components of E independent")
p <- p + theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), title=element_text(size=11.5),
        legend.title=element_blank(), legend.text=element_text(size = 12))

q <- ggplot(resTot.NULL, aes(x=n, y=s.active, fill=meth))+geom_violin(scale = "width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("n") + ylab("Size of active set") + ggtitle("Size of estimated active set of f with p=300, q=5, s=4, components of E independent")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

grid.arrange(p, q, nrow=2)

# E has Toeplitz correlation with rho = 0.4
resTrim.rho04 <- data.frame(matrix(ncol = 4, nrow = nrep * length(n.vec)))
colnames(resTrim.rho04) <- c("n", "MSE", "s.active", "meth")

for(j in 1: length(n.vec)){
  for(i in 1:nrep){
    resTrim.rho04$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim.rho04$n[(j-1)*nrep + i] <- n.vec[j]
    resTrim.rho04$MSE[(j-1)*nrep + i] <- l.VarN.trim.rho04[[j]][[i]][1]
    resTrim.rho04$s.active[(j-1)*nrep + i] <- l.VarN.trim.rho04[[j]][[i]][2]
  }
}

resNone.rho04 <- data.frame(matrix(ncol = 4, nrow = nrep * length(n.vec)))
colnames(resNone.rho04) <- c("n", "MSE", "s.active", "meth")

for(j in 1: length(n.vec)){
  for(i in 1:nrep){
    resNone.rho04$meth[(j-1)*nrep + i] <- "naive"
    resNone.rho04$n[(j-1)*nrep + i] <- n.vec[j]
    resNone.rho04$MSE[(j-1)*nrep + i] <- l.VarN.none.rho04[[j]][[i]][1]
    resNone.rho04$s.active[(j-1)*nrep + i] <- l.VarN.none.rho04[[j]][[i]][2]
  }
}

resTrim.rho04$n <- factor(resTrim.rho04$n, levels = n.vec)
resTrim.rho04$meth <- as.factor(resTrim.rho04$meth)
resNone.rho04$n <- factor(resNone.rho04$n, levels = n.vec)
resNone.rho04$meth <- as.factor(resNone.rho04$meth)

resTot.rho04 <- rbind(resTrim.rho04, resNone.rho04)

p <- ggplot(resTot.rho04, aes(x=n, y=MSE, fill=meth))+geom_violin(scale = "width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("n") + ylab("MSE") + ggtitle("MSE of f with p=300, q=5, s=4, E Toeplitz(0.4)")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

q <- ggplot(resTot.rho04, aes(x=n, y=s.active, fill=meth))+geom_violin(scale = "width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("n") + ylab("Size of active set") + ggtitle("Size of estimated active set of f with p=300, q=5, s=4, E Toeplitz(0.4)")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

grid.arrange(p, q, nrow=2)

# E has Toeplitz correlation with rho = 0.8
resTrim.rho08 <- data.frame(matrix(ncol = 4, nrow = nrep * length(n.vec)))
colnames(resTrim.rho08) <- c("n", "MSE", "s.active", "meth")

for(j in 1: length(n.vec)){
  for(i in 1:nrep){
    resTrim.rho08$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim.rho08$n[(j-1)*nrep + i] <- n.vec[j]
    resTrim.rho08$MSE[(j-1)*nrep + i] <- l.VarN.trim.rho08[[j]][[i]][1]
    resTrim.rho08$s.active[(j-1)*nrep + i] <- l.VarN.trim.rho08[[j]][[i]][2]
  }
}

resNone.rho08 <- data.frame(matrix(ncol = 4, nrow = nrep * length(n.vec)))
colnames(resNone.rho08) <- c("n", "MSE", "s.active", "meth")

for(j in 1: length(n.vec)){
  for(i in 1:nrep){
    resNone.rho08$meth[(j-1)*nrep + i] <- "naive"
    resNone.rho08$n[(j-1)*nrep + i] <- n.vec[j]
    resNone.rho08$MSE[(j-1)*nrep + i] <- l.VarN.none.rho08[[j]][[i]][1]
    resNone.rho08$s.active[(j-1)*nrep + i] <- l.VarN.none.rho08[[j]][[i]][2]
  }
}

resTrim.rho08$n <- factor(resTrim.rho08$n, levels = n.vec)
resTrim.rho08$meth <- as.factor(resTrim.rho08$meth)
resNone.rho08$n <- factor(resNone.rho08$n, levels = n.vec)
resNone.rho08$meth <- as.factor(resNone.rho08$meth)

resTot.rho08 <- rbind(resTrim.rho08, resNone.rho08)



p <- ggplot(resTot.rho08, aes(x=n, y=MSE, fill=meth))+geom_violin(scale = "width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("n") + ylab("MSE") + ggtitle("MSE of f with p=300, q=5, s=4, E Toeplitz(0.8)")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

q <- ggplot(resTot.rho08, aes(x=n, y=s.active, fill=meth))+geom_violin(scale = "width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("n") + ylab("Size of active set") + ggtitle("Size of estimated active set of f with p=300, q=5, s=4, E Toeplitz(0.8)")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

grid.arrange(p, q, nrow=2)