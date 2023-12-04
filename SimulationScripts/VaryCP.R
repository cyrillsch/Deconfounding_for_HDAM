## Code to reproduce Figure 7
# Vary the proportion of confounded covariates

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


# fix n = 400, p = 500

q <- 5
n <-  400
p <- 500

one.sim <- function(prop, meth, seed.val){
  set.seed(seed.val)
  # generate data
  Gamma <- matrix(runif(q*p, min = -1, max = 1), nrow=q)
  # only prop fraction of entries of Gamma should be nonzero
  Gamma0 <- matrix(rbinom(q*p, 1, prop), nrow=q)
  Gamma <- Gamma*Gamma0
  psi <- runif(q, min = 0, max = 2)
  H <- matrix(rnorm(n*q), nrow = n)
  E <- matrix(rnorm(n*p), nrow = n)
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

nrep <- 50


ta <- Sys.time()
prop.vec <- seq(0, 1, 0.05)

l.trim <- list()
l.none <- list()

set.seed(1432)
smax <- 2100000000
seed.vec <- sample(1:smax, nrep)

for (i in 1:length(prop.vec)){
  prop <- prop.vec[i]
  fun = function(seed.val){return(one.sim(prop, "trim", seed.val))}
  l.trim[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  #l.trim[[i]] <- lapply(1:nrep, fun)
  fun = function(seed.val){return(one.sim(prop, "none", seed.val))}
  l.none[[i]] <- mclapply(seed.vec, fun, mc.cores = 60)
  #l.none[[i]] <- lapply(1:nrep, fun)
}
te <- Sys.time()
time.needed <- te-ta

save(l.trim, l.none, time.needed, file = "ResultsVarConfProp_23_10_11.RData")


load("ResultsVarConfProp_23_10_11.RData")


# Plot Results

resTrim <- data.frame(matrix(ncol = 4, nrow = nrep * length(prop.vec)))
colnames(resTrim) <- c("prop", "MSE", "s.active", "meth")

for(j in 1: length(prop.vec)){
  for(i in 1:nrep){
    resTrim$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim$prop[(j-1)*nrep + i] <- prop.vec[j]
    resTrim$MSE[(j-1)*nrep + i] <- l.trim[[j]][[i]][1]
    resTrim$s.active[(j-1)*nrep + i] <- l.trim[[j]][[i]][2]
  }
}


resNone <- data.frame(matrix(ncol = 4, nrow = nrep * length(prop.vec)))
colnames(resNone) <- c("prop", "MSE", "s.active", "meth")

for(j in 1: length(prop.vec)){
  for(i in 1:nrep){
    resNone$meth[(j-1)*nrep + i] <- "naive"
    resNone$prop[(j-1)*nrep + i] <- prop.vec[j]
    resNone$MSE[(j-1)*nrep + i] <- l.none[[j]][[i]][1]
    resNone$s.active[(j-1)*nrep + i] <- l.none[[j]][[i]][2]
  }
}

resTrim$prop <- factor(resTrim$prop, levels = prop.vec)
resTrim$meth <- as.factor(resTrim$meth)
resNone$prop <- factor(resNone$prop, levels = prop.vec)
resNone$meth <- as.factor(resNone$meth)

resTot <- rbind(resTrim, resNone)


library(gridExtra)
p <- ggplot(resTot, aes(x=prop, y=MSE, fill=meth))+geom_violin(scale="width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("prop") + ylab("MSE") + ggtitle("MSE of f with n=400, p=500, q=5, s=4")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

q <- ggplot(resTot, aes(x=prop, y=s.active, fill=meth))+geom_violin(scale="width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("prop") + ylab("Size of active set") + ggtitle("Size of estimated active set of f with n=400, p=500, q=5, s=4")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

grid.arrange(p, q, nrow=2)
