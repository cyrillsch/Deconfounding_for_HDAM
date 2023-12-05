## Code to reproduce Figures 9, 10, 11, 12
# Analysis of the motif regression data set

# Import the functions needed
source("../FunctionsHDAM/FitDeconfoundedHDAM.R")
source("../FunctionsHDAM/AnalyzeFittedHDAM.R")

# The original motif data set from "M. A. Beer and S. Tavazoie (2004). 
# Predicting Gene Expression from Sequence, Cell, Volume 117, Issue 2"
# can be obtained at 
# https://www.sciencedirect.com/science/article/pii/S0092867404003046?via%3Dihub#aep-section-id22
# However, we analyzed a already pre-processed version of the data set from
# "Z. Guo, W. Yuan and C. Zhang (2019). Local Inference in Additive Models with 
# Decorrelated Local Linear Estimator. arXiv preprint arXiv:1907.12732."
# which is not publicly available at the moment.

load("MotifData/MotifData.RData")

# Data as in DLL-Real-Data.R,(Guo, Z., Yuan W. and Zhang, C. (2019).
# Local Inference in Additive Models with Decorrelated Local Linear Estimator.)
# We take the same response as there, i.e. y = genes[, 131].

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

save(fit.null, fit.trim, file = "MotifEvaluation_2023_10_13.RData")

load("MotifEvaluation_2023_10_13.RData")

length(fit.null$active)
# fit.null has 211 active variables

length(fit.trim$active)
# fit.trim has 95 active variables

lengthsq <- function(v){sum(v^2)}
coeflength.null <- as.numeric(lapply(fit.null$coefs, lengthsq))
coeflength.trim <- as.numeric(lapply(fit.trim$coefs, lengthsq))

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

# Plot the strongest functions (from the naive method) that are set to zero by using the deconfounded method
ord <- order(coeflength.null)
# active set of "none" minus active set of "trim"
active.diff <- setdiff(fit.null$active, fit.trim$active)

foo <- function(j){return(j %in% active.diff)}

in.diff <- sapply(ord, foo)
in.diff

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


# coefficient length of deconfounded vs naive
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