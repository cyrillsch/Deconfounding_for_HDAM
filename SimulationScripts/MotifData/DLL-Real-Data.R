### read the data
load("~/Dropbox/Wei-Yuan/Decorrelated Local Linear Estimator/Motif-Regression-Larger/fornicolai.RData")
library(SIHR)
### outcome variables
# summary(genes)
####### 2587 genes/observations
####### 255 columns corresponding to 255 conditions
####### the first 173 columns have their corresponding names
### motifs
# summary(motifs)
####### 2587 genes/observations
####### 666 columns correpsonding to 666 motifs
####### there is a data file explaining what are the motifs


# data
t = 131
y = genes[,t]
X = motifs


DLL.est = matrix(NA,10,3)
DLL.est.se = matrix(NA,10,3)
LF.mat = matrix(NA,10,5)

E = diag(ncol(X))

# foi = c(1,3,13,16,37,41,53,87,89,439)
foi = 16

for (i in 1:length(foi)) {
  print(i)
  mean.i = mean(X[,foi[i]])
  se.i = sd(X[,foi[i]])
  x.eval = c(mean.i,mean.i-se.i,mean.i+se.i)
  DLL.model = DLL(X,y,x.eval=x.eval,foi.ind = foi[i],quant.trans = FALSE,data.swap=FALSE,bwmethod = "rot")
  DLL.est[i,] = DLL.model$est
  DLL.est.se[i,] = DLL.model$est.se
  
  LF.model = LF(X = X, y = y, loading = E[foi[i],])
  LF.mat[i,1] = LF.model$prop.est
  LF.mat[i,2] = LF.model$se
  LF.mat[i,3:4] = LF.model$CI
  LF.mat[i,5] = LF.model$sd.est
  
}


DLL.CI.l = DLL.est - 1.96*DLL.est.se
DLL.CI.u = DLL.est + 1.96*DLL.est.se


mean(LF.mat[,4]-LF.mat[,3])
3.92*apply(DLL.est.se,2,mean)
mean(LF.mat[,5])


