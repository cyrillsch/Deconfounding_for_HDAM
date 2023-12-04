## Code to reproduce Figure 1
# Introductory example

# The R-Script VarP.R needs to be run before.

# consider the setting n = 300, p = 800 from the simulation script VarP.R
load("VarP_2023_10_12.RData")

nrep <- 100
p.vec <- c(50, 100, 200, 400, 800)

resTrim.NULL <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resTrim.NULL) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resTrim.NULL$meth[(j-1)*nrep + i] <- "trim"
    resTrim.NULL$p[(j-1)*nrep + i] <- p.vec[j]
    resTrim.NULL$MSE[(j-1)*nrep + i] <- l.VarP.trim.NULL[[j]][[i]][1]
    resTrim.NULL$s.active[(j-1)*nrep + i] <- l.VarP.trim.NULL[[j]][[i]][2]
  }
}

resNone.NULL <- data.frame(matrix(ncol = 4, nrow = nrep * length(p.vec)))
colnames(resNone.NULL) <- c("p", "MSE", "s.active", "meth")

for(j in 1: length(p.vec)){
  for(i in 1:nrep){
    resNone.NULL$meth[(j-1)*nrep + i] <- "none"
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

# plot histograms
resTrim.p800 <- resTrim.NULL[which(resTrim.NULL$p == "800"),]
resNone.p800 <- resNone.NULL[which(resNone.NULL$p == "800"),]

par(mfrow = c(1,2))
hist(resTrim.p800$MSE, xlim = c(0, 12.5), breaks = seq(0, 12.5, 0.25),
     col = rgb(0,0, 1,0.2), xlab = "MSE", main = "Prediction error of f")
hist(resNone.p800$MSE, add = TRUE, breaks = seq(0, 12.5, 0.25), col = rgb(1,0, 0,0.2))

legend("topright", legend=c("deconfounded","naive"), col=c(rgb(0,0,1,0.2), 
                                                      rgb(1,0,0,0.2)), pt.cex=2, pch=15 )

hist(resTrim.p800$s.active, xlim = c(0, 145), ylim = c(0, 22), breaks = seq(0, 145, 5),
     col = rgb(0,0, 1,0.2), xlab = "Size of estimated active set", main = "Size of estimated active set")
hist(resNone.p800$s.active, add = TRUE, breaks = seq(0, 145, 5), col = rgb(1,0, 0,0.2))
legend("topright", legend=c("deconfounded","naive"), col=c(rgb(0,0,1,0.2), 
                                                   rgb(1,0,0,0.2)), pt.cex=2, pch=15 )
abline(v=4)