source("Imports.R")
library("parallel")
library(ggplot2)
detectCores()

f1 <- function(x){-sin(2*x)}
f2 <- function(x){2-2*tanh(x+0.5)}
f3 <- function(x){x}
f4 <- function(x){4/(exp(x)+exp(-x))}
f <- function(X){f1(X[,1])+f2(X[,2])+f3(X[,3])+f4(X[,4])}

# fix q = 5

q <- 5


one.sim <- function(n, p, meth,  seed.val){
  # generate data
  set.seed(seed.val)
  Gamma <- matrix(runif(q*p, min = -1, max = 1), nrow=q)
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

nrep <- 25



n.vec <- c(250, 500, 750, 1000, 1250, 1500)
p.vec <- c(250, 500, 750, 1000, 1250, 1500)


l.trim <- list()
l.none <- list()


set.seed(1432)
smax <- 2100000000
seed.vec <- sample(1:smax, nrep)

for (j in 1:length(p.vec)){
  for (i in 1:length(n.vec)){
    ta <- Sys.time()
    n <- n.vec[[i]]
    p <- p.vec[[j]]
    fun = function(seed.val){return(one.sim(n, p, "trim", seed.val = seed.val))}
    # this is superfluent. We would not need l.trim[[i]]. l.trim would be enough
    l.trim[[i]] <- mclapply(seed.vec, fun, mc.cores = 20)
    fun = function(seed.val){return(one.sim(n,p, "none", seed.val = seed.val))}
    # this is superfluent. We would not need l.none[[i]]. l.none would be enough
    l.none[[i]] <- mclapply(seed.vec, fun, mc.cores = 20)
    te <- Sys.time()
    time.needed <- te-ta
    filename <- paste("VarNPResults/R_2023_10_27_N",n,"P", p,".RData", sep= "")
    save(l.trim, l.none, time.needed, file = filename)
  }
 
}

# evaluate runtime
rt <- matrix(NA, nrow = length(n.vec), ncol = length(p.vec))
for(i in 1:length(n.vec)){
  for(j in 1:length(p.vec)){
    n <- n.vec[[i]]
    p <- p.vec[[j]]
    filename <- paste("VarNPResults/R_2023_10_27_N",n,"P", p,".RData", sep= "")
    load(filename)
    rt[i,j] <- time.needed
  }
}

rt
# the first few entries are in minutes rather than hours
for(j in 1:5){
  rt[1, j] <- rt[1,j]/60
}
rt[2,1] <- rt[2,1]/60
rt

# dependence on n
par(mfrow=c(2,3))
  for(j in 1:6){
    plot(n.vec, rt[,j])
  }

# dependence on p
par(mfrow=c(2,3))
for(i in 1:6){
  plot(p.vec, rt[i,])
}

# Matrix with mean MSEs and active size
# evaluate runtime
mMSE.trim <- matrix(NA, nrow = length(n.vec), ncol = length(p.vec))
mMSE.none <- matrix(NA, nrow = length(n.vec), ncol = length(p.vec))
mAct.trim <- matrix(NA, nrow = length(n.vec), ncol = length(p.vec))
mAct.none <- matrix(NA, nrow = length(n.vec), ncol = length(p.vec))
for(i in 1:length(n.vec)){
  for(j in 1:length(p.vec)){
    n <- n.vec[[i]]
    p <- p.vec[[j]]
    filename <- paste("VarNPResults/R_2023_10_27_N",n,"P", p,".RData", sep= "")
    load(filename)
    # by typo, we have 6 lists instead of 1. We need to only consider the ith list.
    res.trim <- matrix(unlist(l.trim[[i]]), ncol = 2, byrow = TRUE)
    res.none <- matrix(unlist(l.none[[i]]), ncol = 2, byrow = TRUE)
    mMSE.trim[i,j] <- apply(res.trim, 2, mean)[1]
    mMSE.none[i,j] <- apply(res.none, 2, mean)[1]
    mAct.trim[i,j] <- apply(res.trim, 2, mean)[2]
    mAct.none[i,j] <- apply(res.none, 2, mean)[2]
  }
}

# dependence on n
par(mfrow=c(2,3))
for(j in 1:6){
  plot(n.vec, mMSE.trim[, j])
}

# dependence on p
par(mfrow=c(2,3))
for(i in 1:6){
  plot(p.vec, mMSE.trim[i,])
}


# generate "long" dataframe to make linear regression.
MSE.long <- numeric(length(n.vec)*length(p.vec))
n.long <- rep(n.vec, length(p.vec))
p.long <- rep(p.vec, each = length(n.vec))

for(j in 1:length(p.vec)){
  MSE.long[((j-1)*length(n.vec)+1):(j*length(n.vec))] <- mMSE.trim[,j]
}

dat <- data.frame(MSE=MSE.long, n = n.long, p= p.long)

par(mfrow=c(1,1))
plot(dat$n, dat$MSE)
plot(log(dat$n), log(dat$MSE))

fit <- lm(log(MSE)~log(n), data = dat)
summary(fit)
abline(fit)

# The convergence seems to be faster than it should actually be
