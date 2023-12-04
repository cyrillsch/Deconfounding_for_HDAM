## Code to reproduce Figures 4 and 5
# Vary the dimension p of X

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

# fix n = 300, q = 5

q <- 5
n <- 300

one.sim <- function(p, meth, rho = NULL, seed.val){
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
p.vec <- c(50, 100, 200, 400, 800)


l.VarP.trim.NULL <- list()
l.VarP.none.NULL <- list()
l.VarP.trim.rho04 <- list()
l.VarP.none.rho04 <- list()
l.VarP.trim.rho08 <- list()
l.VarP.none.rho08 <- list()

set.seed(1432)
smax <- 2100000000
seed.vec <- sample(1:smax, nrep)

for (i in 1:length(p.vec)){
  p <- p.vec[[i]]
  rho = NULL
  fun = function(seed.val){return(one.sim(p, "trim", rho = NULL, seed.val = seed.val))}
  l.VarP.trim.NULL[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  fun = function(seed.val){return(one.sim(p, "none", rho = NULL, seed.val = seed.val))}
  l.VarP.none.NULL[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  # rho = 0.4
  fun = function(seed.val){return(one.sim(p, "trim", rho = 0.4, seed.val = seed.val))}
  l.VarP.trim.rho04[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  fun = function(seed.val){return(one.sim(p, "none", rho = 0.4, seed.val = seed.val))}
  l.VarP.none.rho04[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  # rho = 0.8
  fun = function(seed.val){return(one.sim(p, "trim", rho = 0.8, seed.val = seed.val))}
  l.VarP.trim.rho08[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  fun = function(seed.val){return(one.sim(p, "none", rho = 0.8, seed.val = seed.val))}
  l.VarP.none.rho08[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
}
te <- Sys.time()
time.needed <- te-ta

 save(l.VarP.trim.NULL, l.VarP.none.NULL, l.VarP.trim.rho04, l.VarP.none.rho04, l.VarP.trim.rho08, l.VarP.none.rho08, time.needed, file = "VarP_2023_10_12.RData")

# Plot results
load("VarP_2023_10_12.RData")


# Without correlation in E
resTrim.NULL <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resTrim.NULL) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resTrim.NULL$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim.NULL$p[(j-1)*nrep + i] <- p.vec[j]
    resTrim.NULL$MSE[(j-1)*nrep + i] <- l.VarP.trim.NULL[[j]][[i]][1]
    resTrim.NULL$s.active[(j-1)*nrep + i] <- l.VarP.trim.NULL[[j]][[i]][2]
  }
}

resNone.NULL <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resNone.NULL) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resNone.NULL$meth[(j-1)*nrep + i] <- "naive"
    resNone.NULL$p[(j-1)*nrep + i] <- p.vec[j]
    resNone.NULL$MSE[(j-1)*nrep + i] <- l.VarP.none.NULL[[j]][[i]][1]
    resNone.NULL$s.active[(j-1)*nrep + i] <- l.VarP.none.NULL[[j]][[i]][2]
  }
}

resTrim.NULL$p <- factor(resTrim.NULL$p, levels = p.vec)
resTrim.NULL$meth <- as.factor(resTrim.NULL$meth)
resNone.NULL$p <- factor(resNone.NULL$p, levels = p.vec)
resNone.NULL$meth <- as.factor(resNone.NULL$meth)

resTot.NULL <- rbind(resTrim.NULL, resNone.NULL)



library(gridExtra)
p <- ggplot(resTot.NULL, aes(x=p, y=MSE, fill=meth))+geom_violin(scale = "width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("p") + ylab("MSE") + ggtitle("MSE of f with n=300, q=5, s=4, components of E independent")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

q <- ggplot(resTot.NULL, aes(x=p, y=s.active, fill=meth))+geom_violin(scale = "width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("p") + ylab("Size of active set") + ggtitle("Size of estimated active set of f with n=300, q=5, s=4, components of E independent")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

grid.arrange(p, q, nrow=2)



# E has Toeplitz correlation with rho = 0.4
resTrim.rho04 <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resTrim.rho04) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resTrim.rho04$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim.rho04$p[(j-1)*nrep + i] <- p.vec[j]
    resTrim.rho04$MSE[(j-1)*nrep + i] <- l.VarP.trim.rho04[[j]][[i]][1]
    resTrim.rho04$s.active[(j-1)*nrep + i] <- l.VarP.trim.rho04[[j]][[i]][2]
  }
}

resNone.rho04 <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resNone.rho04) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resNone.rho04$meth[(j-1)*nrep + i] <- "naive"
    resNone.rho04$p[(j-1)*nrep + i] <- p.vec[j]
    resNone.rho04$MSE[(j-1)*nrep + i] <- l.VarP.none.rho04[[j]][[i]][1]
    resNone.rho04$s.active[(j-1)*nrep + i] <- l.VarP.none.rho04[[j]][[i]][2]
  }
}

resTrim.rho04$p <- factor(resTrim.rho04$p, levels = p.vec)
resTrim.rho04$meth <- as.factor(resTrim.rho04$meth)
resNone.rho04$p <- factor(resNone.rho04$p, levels = p.vec)
resNone.rho04$meth <- as.factor(resNone.rho04$meth)

resTot.rho04 <- rbind(resTrim.rho04, resNone.rho04)



p <- ggplot(resTot.rho04, aes(x=p, y=MSE, fill=meth))+geom_violin(scale = "width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("p") + ylab("MSE") + ggtitle("MSE of f with n=300, q=5, s=4, E Toeplitz(0.4)")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

q <- ggplot(resTot.rho04, aes(x=p, y=s.active, fill=meth))+geom_violin(scale = "width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("p") + ylab("Size of active set") + ggtitle("Size of estimated active set of f with n=300, q=5, s=4, E Toeplitz(0.4)")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))
grid.arrange(p, q, nrow=2)


# E has Toeplitz correlation with rho = 0.8
resTrim.rho08 <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resTrim.rho08) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resTrim.rho08$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim.rho08$p[(j-1)*nrep + i] <- p.vec[j]
    resTrim.rho08$MSE[(j-1)*nrep + i] <- l.VarP.trim.rho08[[j]][[i]][1]
    resTrim.rho08$s.active[(j-1)*nrep + i] <- l.VarP.trim.rho08[[j]][[i]][2]
  }
}

resNone.rho08 <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resNone.rho08) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resNone.rho08$meth[(j-1)*nrep + i] <- "naive"
    resNone.rho08$p[(j-1)*nrep + i] <- p.vec[j]
    resNone.rho08$MSE[(j-1)*nrep + i] <- l.VarP.none.rho08[[j]][[i]][1]
    resNone.rho08$s.active[(j-1)*nrep + i] <- l.VarP.none.rho08[[j]][[i]][2]
  }
}

resTrim.rho08$p <- factor(resTrim.rho08$p, levels = p.vec)
resTrim.rho08$meth <- as.factor(resTrim.rho08$meth)
resNone.rho08$p <- factor(resNone.rho08$p, levels = p.vec)
resNone.rho08$meth <- as.factor(resNone.rho08$meth)

resTot.rho08 <- rbind(resTrim.rho08, resNone.rho08)



p <- ggplot(resTot.rho08, aes(x=p, y=MSE, fill=meth))+geom_violin(scale = "width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("p") + ylab("MSE") + ggtitle("MSE of f with n=300, q=5, s=4, E Toeplitz(0.8)")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

q <- ggplot(resTot.rho08, aes(x=p, y=s.active, fill=meth))+geom_violin(scale = "width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("p") + ylab("Size of active set") + ggtitle("Size of estimated active set of f with n=300, q=5, s=4, E Toeplitz(0.8)")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))
grid.arrange(p, q, nrow=2)



