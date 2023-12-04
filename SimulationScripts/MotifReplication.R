# Check, how reasonable the results from the motiv example are
# We extract the first few principal components of the data and construct
# a synthetic response using them.

source("Imports.R")
load("MotifData/fornicolai.RData")


#Data as in DLL-Real-Data.R, i.e. y=genes[, 131].
# data
t = 131
y = genes[,t]
X = motifs
par(mfrow=c(1,1))
svdX.d <- svd(X)$d
plot(svdX.d)
# Very large first eigenvalue. But Data not centered.
Xmeans <- colMeans(X)
Xc <- X-outer(rep(1, length(y)), Xmeans)
svdX.c <- svd(Xc)
plot(svdX.c$d, main = "Singular values of (centered) motif data", ylab="singular value d_l", xlab = "l")
# still confounding might be present
abline(v = 20)




# take the first 20 principal components as the confounder H
n <- length(y)
H <- sqrt(n)*svdX.c$u[, 1:20]

# generate synthetic Y

f <- function(a, t){0.8*sin(1.5*(t-a))*exp(0.8*t)}

# With a  randomly, this yields functions similar to the
# ones fitted on the Motif data.

set.seed(1220)
par(mfrow=c(3, 3))
xx <- seq(-1.87, 1.54, 0.01)
for(i in 1:9){
  a <- runif(1)
  plot(xx, f(a, xx), type = "l")
}


# sample 95 active variables at random

act <- sample(1:666, 95, replace = FALSE)

# Only the first few functions should be "strong". I.e. weights decreasing and random sign
Y <- 0
# w <-  1/sqrt(1:length(act))*(2*rbinom(length(act), 1, 0.5)-1)
w <- -c(1, -1, -1, 1, 1/1:(length(act)-4)*(2*rbinom(length(act)-4, 1, 0.5)-1))
a <- runif(length(act))
for(i in 1:length(act)){
  Y <- Y+w[i]*f(a[i], X[, act[i]])
}
ftrue <- Y

# # without confounding in Y
# Y <- Y+0.5*rnorm(n)

# with confounding in Y
del <- 0.1*runif(20, -1,1)
Y <- Y + c(H%*%del)+0.5*rnorm(n)

indi <- sample(1:666, 9)
par(mfrow=c(3, 3))
for(i in 1:9){
  plot(Xc[, indi[i]], Y)
}

for(i in 1:9){
  plot(Xc[, indi[i]], y)
}



fit.null <- FitDeconfoundedHDAM(Y, Xc, n.K = 5, meth = "none", cv.method = "1se", cv.k = 5, n.lambda1 = 15, n.lambda2 = 50)
fit.trim <- FitDeconfoundedHDAM(Y, Xc, n.K = 5, meth = "trim", cv.method = "1se", cv.k = 5, n.lambda1 = 15, n.lambda2 = 50)

#save(fit.null, fit.trim, file = "MotifReplicationWithConfounding_2023_10_23.RData")

load("MotifReplicationWithConfounding_2023_10_23.RData")

length(fit.null$active)
# fit.null has 100 active variables

length(fit.trim$active)
# fit.trim has 53 active variables

intersect(act, fit.trim$active)
intersect(act, fit.null$active)
act
# both, fit.trim and fit.null capture the most important active variables from act.

lengthsq <- function(v){sum(v^2)}
coeflength.null <- as.numeric(lapply(fit.null$coefs, lengthsq))
coeflength.trim <- as.numeric(lapply(fit.trim$coefs, lengthsq))

par(mfrow=c(2,1))
plot(coeflength.null, ylim = c(0, 10), main = "Length of coefficient vector without transformation", xlab = "j", ylab = "length squared")
plot(coeflength.trim, ylim = c(0,10), main = "Length of coefficient vector witwith trim transformation", xlab = "j", ylab = "length squared")

par(mfrow=c(1,1))
plot(coeflength.null, coeflength.trim)

# the coeficient lengths are very similar.


#plot 9 strongest functions in terms of coeflength

ord <- order(coeflength.trim)
ind <- tail(ord, n=9)



par(mfrow = c(3,3))
for (l in 1:9){
  j <- ind[10-l]
  xx <- seq(min(Xc[,j]), max(Xc[,j]), length.out=50)
  fhatj.trim <- estimate.fj(xx, j, fit.trim)
  plot(xx, fhatj.trim, type = "l", col = "blue", xlim=c(min(xx),max(xx)), xlab = "X_j", ylab = "f_j(X_j)", main = paste("j=", j), ylim = c(-0.5, 2.5))
  fhatj.null <- estimate.fj(xx, j, fit.null)
  lines(xx, fhatj.null, type = "l", col = "red", ylim=c(min(xx),max(xx)), xlab = "X_j", ylab = "f_j(X_j)", main = paste("j=", j))
  legend("topleft", legend = c("trim transform", "no transformation"), col = c("blue", "red"), lwd = 1)

  abline(h=0, col="grey")
  points(Xc[,j], rep(-0.5, length(X[, j])), pch = 16, cex= 1, col=rgb(0,0,0, 0.1))
}

# fitted functions look very similar
# mtext("Fitted functions f_j for selected j", side = 3, line = - 2, cex = 1.5, outer = TRUE)


fhat.null <- estimate.function(Xc, fit.null)
fhat.trim <- estimate.function(Xc, fit.trim)



mean(fhat.null^2)
mean(fhat.trim^2)

# We know the "true" fitted values
mean(ftrue^2)

mean(fhat.null^2)/mean(fhat.trim^2)

# fhat.null and fhat.trim have very similar size. 
# Also very similar to ftrue

par(mfrow=c(1,1))
plot(fhat.null, fhat.trim)
fit <- lm(fhat.trim~fhat.null)
summary(fit)
abline(fit, col = "red")

plot(fhat.trim, fhat.null)
fit <- lm(fhat.null~fhat.trim)
abline(fit, col = "red")
summary(fit)



# The resulted fhat.trim and fhat.null are very similar to the ftrue
# This result is quite different

