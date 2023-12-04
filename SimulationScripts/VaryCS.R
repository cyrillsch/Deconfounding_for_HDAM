# Vary the strength of the confounding

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

one.sim <- function(cs, meth, seed.val){
  set.seed(seed.val)
  # generate data
  Gamma <- matrix(runif(q*p, min = -1, max = 1), nrow=q)
  psi <- runif(q, min = 0, max = cs)
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
set.seed(1527)


ta <- Sys.time()
cs.vec <- seq(0, 3, by=0.25)

l.trim <- list()
l.none <- list()

set.seed(1432)
smax <- 2100000000
seed.vec <- sample(1:smax, nrep)

for (i in 1:length(cs.vec)){
  cs <- cs.vec[i]
  fun = function(seed.val){return(one.sim(cs, "trim", seed.val))}
  l.trim[[i]] <- mclapply(seed.vec, fun, mc.cores = 20)
  #l.trim[[i]] <- lapply(1:nrep, fun)
  fun = function(seed.val){return(one.sim(cs, "none", seed.val))}
  l.none[[i]] <- mclapply(seed.vec, fun, mc.cores = 20)
  #l.none[[i]] <- lapply(1:nrep, fun)
}
te <- Sys.time()
time.needed <- te-ta

# save(l.trim, l.none, time.needed, file = "ResultsVarCS_23_10_11.RData")

load("ResultsVarCS_23_10_11.RData")


resTrim <- data.frame(matrix(ncol = 4, nrow = nrep * length(cs.vec)))
colnames(resTrim) <- c("cs", "MSE", "s.active", "meth")

for(j in 1: length(cs.vec)){
  for(i in 1:nrep){
    resTrim$meth[(j-1)*nrep + i] <- "deconfounded"
    resTrim$cs[(j-1)*nrep + i] <- cs.vec[j]
    resTrim$MSE[(j-1)*nrep + i] <- l.trim[[j]][[i]][1]
    resTrim$s.active[(j-1)*nrep + i] <- l.trim[[j]][[i]][2]
  }
}

resNone <- data.frame(matrix(ncol = 4, nrow = nrep * length(cs.vec)))
colnames(resNone) <- c("cs", "MSE", "s.active", "meth")

for(j in 1: length(cs.vec)){
  for(i in 1:nrep){
    resNone$meth[(j-1)*nrep + i] <- "naive"
    resNone$cs[(j-1)*nrep + i] <- cs.vec[j]
    resNone$MSE[(j-1)*nrep + i] <- l.none[[j]][[i]][1]
    resNone$s.active[(j-1)*nrep + i] <- l.none[[j]][[i]][2]
  }
}

resTrim$cs <- factor(resTrim$cs, levels = cs.vec)
resTrim$meth <- as.factor(resTrim$meth)
resNone$cs <- factor(resNone$cs, levels = cs.vec)
resNone$meth <- as.factor(resNone$meth)

resTot <- rbind(resTrim, resNone)


library(gridExtra)
p <- ggplot(resTot, aes(x=cs, y=MSE, fill=meth))+geom_violin(scale="width")
p <- p + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("confounding strength") + ylab("MSE") + ggtitle("MSE of f with n=400, p=500, q=5, s=4")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))



q <- ggplot(resTot, aes(x=cs, y=s.active, fill=meth))+geom_violin(scale="width")
q <- q + stat_summary(fun.y=mean, geom="point", position=position_dodge(0.9)) +
  xlab("confounding strength") + ylab("Size of the active set") + ggtitle("Size of the estimated active set with n=400, p=500, q=5, s=4")
q <- q + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12), title=element_text(size=11.5),
               legend.title=element_blank(), legend.text=element_text(size = 12))

grid.arrange(p, q, nrow=2)
