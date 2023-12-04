
source("Imports.R")

f1 <- function(x){-sin(2*x)}
f2 <- function(x){2-2*tanh(x+0.5)}
f3 <- function(x){x}
f4 <- function(x){4/(exp(x)+exp(-x))}
f <- function(X){f1(X[,1])+f2(X[,2])+f3(X[,3])+f4(X[,4])}

# fix p = 300, q = 5

q <- 5
p <- 300
n <- 600


# generate data
set.seed(855)
Gamma <- matrix(runif(q*p, min = -1, max = 1), nrow=q)
psi <- runif(q, min = 0, max = 2)
H <- matrix(rnorm(n*q), nrow = n)


E <- matrix(rnorm(n*p), nrow = n)

e <- 0.5*rnorm(n)
X <- H %*% Gamma + E
Y <- f(X) + H %*% psi + e






fit.null <- FitDeconfoundedHDAM(Y, X, n.K = 5, meth = "none", cv.method = "1se", cv.k = 5, n.lambda1 = 15, n.lambda2 = 25)
fit.trim <- FitDeconfoundedHDAM(Y, X, n.K = 5, meth = "trim", cv.method = "1se", cv.k = 5, n.lambda1 = 15, n.lambda2 = 25)

#save(fit.null, fit.trim, file = "MotifReplicationWithConfounding_2023_10_23.RData")



length(fit.null$active)
# fit.null has 100 active variables

length(fit.trim$active)
# fit.trim has 53 active variables


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



par(mfrow = c(2,2))
for (l in 1:9){
  j <- ind[10-l]
  xx <- seq(min(X[,j]), max(X[,j]), length.out=50)
  fhatj.trim <- estimate.fj(xx, j, fit.trim)
  plot(xx, fhatj.trim, type = "l", col = "blue", xlim=c(min(xx),max(xx)), xlab = "X_j", ylab = "f_j(X_j)", main = paste("j=", j), ylim = c(-0.5, 2.5))
  fhatj.null <- estimate.fj(xx, j, fit.null)
  lines(xx, fhatj.null, type = "l", col = "red", ylim=c(min(xx),max(xx)), xlab = "X_j", ylab = "f_j(X_j)", main = paste("j=", j))
  legend("topleft", legend = c("trim transform", "no transformation"), col = c("blue", "red"), lwd = 1)
  
  abline(h=0, col="grey")
  points(X[,j], rep(-0.5, length(X[, j])), pch = 16, cex= 1, col=rgb(0,0,0, 0.1))
}

# fitted functions look very similar
# mtext("Fitted functions f_j for selected j", side = 3, line = - 2, cex = 1.5, outer = TRUE)


fhat.null <- estimate.function(X, fit.null)
fhat.trim <- estimate.function(X, fit.trim)



mean(fhat.null^2)
mean(fhat.trim^2)


# We know the "true" fitted values
ftrue <- f(X)
mean(ftrue^2)

mean(fhat.null^2)/mean(fhat.trim^2)



par(mfrow=c(1,1))
plot(fhat.null, fhat.trim)
fit <- lm(fhat.trim~fhat.null)
summary(fit)
abline(fit, col = "red")

plot(fhat.trim, fhat.null)
fit <- lm(fhat.null~fhat.trim)
abline(fit, col = "red")
summary(fit)


plot(ftrue, fhat.trim)
plot(ftrue, fhat.null)

# fhat.trim fits better, but there is not a non-linear relationship between fhat.trim and fhat.none