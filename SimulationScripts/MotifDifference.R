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

par(mfrow = c(3,3))
for (l in 1:9){
  j <- int.ind[10-l]
  xx <- seq(min(X[,j]), max(X[,j]), length.out=50)
  fhatj.trim <- estimate.fj(xx, j, fit.trim)
  plot(xx, fhatj.trim, type = "l", col = "blue", xlim=c(min(xx),max(xx)), xlab = "X_j", ylab = "f_j(X_j)", main = paste("j=", j), ylim = c(-0.5, 2.5))
  fhatj.null <- estimate.fj(xx, j, fit.null)
  lines(xx, fhatj.null, type = "l", col = "red", ylim=c(min(xx),max(xx)), xlab = "X_j", ylab = "f_j(X_j)", main = paste("j=", j))
  legend("topleft", legend = c("trim transformation", "no transformation"), col = c("blue", "red"), lwd = 1)
  
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
par(mfrow = c(1,2))
plot(sqrt(coeflength.null), sqrt(coeflength.trim), main = "Norm of coefficient vectors", xlab = "no transformation", ylab = "trim transformation")
abline(0,1)
plot(jaccard.sim, main = "Jaccard similarity of top l index sets", type = "l", xlab = "l", ylab = "Jaccard similarity")


