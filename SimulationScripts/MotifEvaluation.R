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
svdX.c <- svd(Xc)$d
plot(svdX.c, main = "Singular values of (centered) motif data", ylab="singular value d_l", xlab = "l")
# still confounding might be present


set.seed(1443)

fit.null <- FitDeconfoundedHDAM(y, X, n.K = 5, meth = "none", cv.method = "1se", cv.k = 5, n.lambda1 = 15, n.lambda2 = 50)
fit.trim <- FitDeconfoundedHDAM(y, X, n.K = 5, meth = "trim", cv.method = "1se", cv.k = 5, n.lambda1 = 15, n.lambda2 = 50)

# save(fit.null, fit.trim, file = "MotifEvaluation_2023_10_13.RData")



load("MotifEvaluation_2023_10_13.RData")

length(fit.null$active)
# fit.null has 211 active variables

length(fit.trim$active)
# fit.trim has 95 active variables

lengthsq <- function(v){sum(v^2)}
coeflength.null <- as.numeric(lapply(fit.null$coefs, lengthsq))
coeflength.trim <- as.numeric(lapply(fit.trim$coefs, lengthsq))

par(mfrow=c(2,1))
plot(coeflength.null, ylim = c(0, 10), main = "Length of coefficient vector without transformation", xlab = "j", ylab = "length squared")
plot(coeflength.trim, ylim = c(0,10), main = "Length of coefficient vector with trim transformation", xlab = "j", ylab = "length squared")

par(mfrow=c(1,1))
plot(coeflength.null, coeflength.trim)
abline(0,1)


#plot 9 strongest functions in terms of coeflength

ord <- order(coeflength.trim)
ind <- tail(ord, n=9)



par(mfrow = c(3,3), mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(2, 0.5, 0))
for (l in 1:9){
  j <- ind[10-l]
  xx <- seq(min(X[,j]), max(X[,j]), length.out=50)
  fhatj.trim <- estimate.fj(xx, j, fit.trim)
  plot(xx, fhatj.trim, type = "l", col = "blue", xlim=c(min(xx),max(xx)), xlab = "Xj", ylab = "fj(Xj)", main = paste("j =", j), ylim = c(-0.5, 2.5))
  fhatj.null <- estimate.fj(xx, j, fit.null)
  lines(xx, fhatj.null, type = "l", col = "red", ylim=c(min(xx),max(xx)), xlab = "Xj", ylab = "fj(Xj)", main = paste("j =", j))
  legend("topleft", legend = c("deconfounded", "naive"), col = c("blue", "red"), lwd = 1)

  abline(h=0, col="grey")
  points(X[,j], rep(-0.5, length(X[, j])), pch = 16, cex= 1, col=rgb(0,0,0, 0.1))
}

# mtext("Fitted functions f_j for selected j", side = 3, line = - 2, cex = 1.5, outer = TRUE)

# fitted functions look similar
# Note that the functions are actually empirically centered, as most of the mass lies in the area, where the functions are (slightly) negative


##############################################################################################
# which are the largest coefficients for "none" that are not in the active set for "trim"?

ord <- order(coeflength.null)
# active set of "none" minus active set of "trim"
active.diff <- setdiff(fit.null$active, fit.trim$active)


foo <- function(j){return(j %in% active.diff)}

in.diff <- sapply(ord, foo)
in.diff

ord[in.diff]

coeflength.null[1]

int.ind <- tail(ord[in.diff], 9)

par(mfrow = c(3,3), mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(2, 0.5, 0))
for (l in 1:9){
  j <- int.ind[10-l]
  xx <- seq(min(X[,j]), max(X[,j]), length.out=50)
  fhatj.trim <- estimate.fj(xx, j, fit.trim)
  plot(xx, fhatj.trim, type = "l", col = "blue", xlim=c(min(xx),max(xx)), xlab = "Xj", ylab = "fj(Xj)", main = paste("j =", j), ylim = c(-0.5, 2.5))
  fhatj.null <- estimate.fj(xx, j, fit.null)
  lines(xx, fhatj.null, type = "l", col = "red", ylim=c(min(xx),max(xx)), xlab = "Xj", ylab = "fj(Xj)", main = paste("j =", j))
  legend("topleft", legend = c("deconfounded", "naive"), col = c("blue", "red"), lwd = 1)
  
  abline(h=0, col="grey")
  points(X[,j], rep(-0.5, length(X[, j])), pch = 16, cex= 1, col=rgb(0,0,0, 0.1))
}

###############################################################################

# are the 92 variables selected by trim more or less the 92 strongest variable for the naive method?

top92.null <- head(order(coeflength.null, decreasing = TRUE), 92)

length(intersect(fit.trim$active, top92.null))
# 68

# coefficient length of trim vs naive
par(mfrow=c(1,1))
plot(sqrt(coeflength.null), sqrt(coeflength.trim), main = "Norm of coefficient vectors", xlab = "no transformation", ylab = "trim transformation")
abline(0,1)


# for l = 1,..., 211, record the size of the intersection of the strongest l variables for trim and the strongest l variables for none
ord.null <- order(coeflength.null, decreasing = TRUE)[1:211]
ord.trim <- order(coeflength.trim, decreasing = TRUE)[1:92]

size.intersect <- numeric(211)
size.union <- numeric(211)

par(mfrow = c(1,1))
for(l in 1:211){
  size.intersect[l] <- length(intersect(head(ord.null, l), head(ord.trim, l)))
  size.union[l] <- length(union(head(ord.null, l), head(ord.trim, l)))
}

jaccard.dis <- (size.union-size.intersect)/size.union
jaccard.sim <- size.intersect/size.union

plot(jaccard.sim, type = "l")

# both plots in one
dev.off()
par(mfrow = c(1,2))
plot(sqrt(coeflength.null), sqrt(coeflength.trim), main = "Norm of coefficient vectors", xlab = "naive", ylab = "deconfounded", xlim = c(0, 3.1), ylim = c(0, 3.1))
abline(0,1)
plot(jaccard.sim, main = "Jaccard similarity of top l index sets", type = "l", xlab = "l", ylab = "Jaccard similarity")






################################################################################
# Additional comparisons (not in the paper)
################################################################################





fhat.null <- estimate.function(X, fit.null)
fhat.trim <- estimate.function(X, fit.trim)



mean(fhat.null^2)
mean(fhat.trim^2)

mean(fhat.null^2)/mean(fhat.trim^2)

# fhat.null seems to overfit more.

par(mfrow=c(1,1))
plot(fhat.null, fhat.trim)
fit <- lm(fhat.trim~fhat.null)
summary(fit)
abline(fit, col = "red")

plot(fhat.trim, fhat.null)
fit <- lm(fhat.null~fhat.trim)
abline(fit, col = "red")
summary(fit)

# there seems to be a slightly nonlinear relationship



plot(fhat.null, y)
fit <- lm(y~fhat.null)
abline(fit)
summary(fit)


plot(fhat.trim, y)
fit1 <- lm(y~fhat.trim)
abline(fit1)
summary(fit1)



# residuals vs. fitted values
plot(y-fhat.null, fhat.null)
plot(y-fhat.trim, fhat.trim)


