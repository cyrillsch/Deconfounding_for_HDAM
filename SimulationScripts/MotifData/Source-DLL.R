### Source-DLL.R
### Functions: Implement Decorrelated Local Linear Estimator as developed
###            in the paper Guo and Zhang (2019), for function derivative
###            estimate in high-dimensional sparse additive models


### Required  packages
library(SAM) # sparse additive models
library(glmnet) # lasso estimator
library(nprobust) # bandwidth selection
library(locpol)
library(np) # bandwidth selection


### cv.SAM
### Function: Cross-validated Sparse Additive Model
###           fitting to select the best number of
###           basis functions in smoothing and the
###           best tuning parameter lambda
### Input: X, continuous, the covariates (n by p matrix)
###        y, continuous, the outcome (n by 1 vector)
###        kfold, integer, the numbers of folds for cross validation,
###               default value is 10
###        degree, vector of integers, a sequence of number of basis
###                functions to be considered, default is 3
###        lam.seq, continuous, a sequence of lambdas to be considered,
###                 usually in decreasing order to speed up, default is
###                 0.5^(seq(0, 10, l=200))
### Output: sam.final, the sparse additive model object as in the SAM
###                  package with the best parameters
###         sigma1.sq, the consistent estimate of the variance of the
###                    noise using the mean squared errors
###         X, input covariates
###         Z.hat, the transformed covariates in a compact support,
###                [0, 1] using quantile transformation in our case
cv.SAM = function(X, y, kfold=5, degree=3, lam.seq=NULL, quant.trans) {
  n = nrow(X); p = ncol(X)
  # quantile transformation of X
  if (quant.trans) {
    Z.hat = apply(X, 2, function(x) ecdf(x)(x))
  } else {
    Z.hat = X
  }

  # default lambda sequence, approximately from 0.001 to 1
  if (is.null(lam.seq)) lam.seq = 0.5^seq(0, 10, l=200)

  # cross validation
  len_d = length(degree)
  len_lam = length(lam.seq)
  MSE = matrix(0, len_d, 3)
  colnames(MSE) = c("degree", "lambda", "MSE")
  # break the index into k folds
  folds = cut(seq(1,n),breaks=kfold,labels=FALSE)
  for (i in 1:len_d) {
    mse.lam = rep(0, len_lam)
    for (fold in 1:kfold) {
      testInd = which(folds == fold, arr.ind = TRUE)
      # create train and test data
      Z.test = Z.hat[testInd, ]
      y.test = y[testInd]
      Z.train = Z.hat[-testInd, ]
      y.train = y[-testInd]
      # provide a sequence of lambda to speed up
      sam.fit = samQL(Z.train, y.train, p = degree[i], lambda = lam.seq)
      y.pred = predict(sam.fit, newdata = Z.test)$values
      mse.lam = mse.lam + apply(y.test-y.pred, 2, function(x) mean(x^2))/kfold
    }

    temp.lam = lam.seq[which.min(mse.lam)]
    MSE[i, "degree"] = degree[i]
    MSE[i, "lambda"] = temp.lam
    MSE[i, "MSE"] = min(mse.lam)
  }

  best.ind = which.min(MSE[, "MSE"])
  best.degree = MSE[best.ind, "degree"]
  best.lam = MSE[best.ind, "lambda"]
  sam.final = samQL(Z.hat, y, p = best.degree, lambda = best.lam)
  sigma1.sq = mean((y-predict(sam.final, newdata = Z.hat)$values)^2)

  returnList = list(sam.final = sam.final,
                     sigma1.sq = sigma1.sq,
                     X = X,
                     Z.hat = Z.hat)
  returnList
}


### predict.SAM
### Function: Calculate the estimate of functions given the
###           fitted sparse additive models and the test data
### Input: sam.obj, fitted sparse additive model in cv.SAM
###        Xt, continuous, the test covariates (nt by p matrix)
### Output: f.hat, estimated functions including an intercept
###                (nt by p+1 matrix)
predict.SAM = function(sam.obj, Xt, quant.trans) {
  nt = nrow(Xt)
  p = ncol(Xt)
  X = sam.obj$X
  sam.final = sam.obj$sam.final
  # scale the X using quantile transformation or not
  if (quant.trans) {
    emp.cdf = apply(X, 2, function(x) ecdf(x))
    X.qt = matrix(0, nt, p)
    for (j in 1:p) {
      X.qt[, j] = emp.cdf[[j]](Xt[, j])
    }
  } else {
    X.qt = Xt
  }
  # max-min transformation default in the SAM package
  X.min.rep = matrix(rep(sam.final$X.min,nt),nrow=nt,byrow=T)
  X.ran.rep = matrix(rep(sam.final$X.ran,nt),nrow=nt,byrow=T)
  
  X.scale = (X.qt-X.min.rep)/X.ran.rep
  X.scale = pmax(X.scale,0)
  X.scale = pmin(X.scale,1)


  # calculate each estimated function
  f.hat = matrix(0, nt, p)
  for (j in 1:p) {
    bspline = ns(X.scale[, j], df = sam.final$p, knots = sam.final$knots[, j],
                  Boundary.knots = sam.final$Boundary.knots[, j]) # B-spline matrix for prediction
    ind.j = (j-1)*sam.final$p + c(1:sam.final$p) # index of the target component function coefficients
    f.hat[, j] = bspline %*% sam.final$w[ind.j]
  }
  # add the intercept
  f.hat = cbind(rep(sam.final$intercept, nt), f.hat)
  colnames(f.hat) = c("Intercept", paste("X", 1:p, sep = ""))

  return(f.hat)
}


### LASSO
### Function: computes the LASSO estimator
### -If lambda is given, use glmnet and standard Lasso
### -If lambda is set to the character string "CV", then glmnet with
### lambda selected by cross-validation is used
### -If lambda is not given or is set to NULL, use square root Lasso
Lasso = function(X, y, lambda = NULL, intercept = TRUE) {
  p = ncol(X)
  n = nrow(X)

  htheta = if (is.null(lambda)) {
    lambda = sqrt(qnorm(1 - (0.1 / p)) / n)
    outLas = slim(X, y, lambda = lambda, method = "lq", q = 2,
                   verbose = FALSE)
    # Objective : sqrt(RSS/n) + lambda * penalty
    c(as.vector(outLas$intercept), as.vector(outLas$beta))
  } else if (lambda == "CV") {
    outLas = cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas = cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "scalreg") {
    Xc = if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas = scalreg(Xc, y)
    # return object
    if (intercept) {
      outLas$coefficients
    } else {
      # add a coefficient for the (not estimated) intercept b/c of implementation
      c(0, outLas$coefficients)
    }
  } else {
    outLas = glmnet(X, y, family = "gaussian", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }

  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}


### Initialization.step
### Function: Computes the initial LASSO estimator and quantities based thereon
### lambda is set as "CV"
Initialization.step = function(X, y, lambda = NULL, intercept = FALSE) {
  n = nrow(X)
  # col.norm = 1 / sqrt((1 / n) * diag(t(X) %*% X))
  col.norm = 1 / sqrt((1 / n) * diag(t(X)%*%X))
  Xnor = X %*% diag(col.norm)

  ### Call Lasso
  htheta = Lasso(Xnor, y, lambda = lambda, intercept = intercept)

  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb = cbind(rep(1, n), Xnor)
    col.norm = c(1, col.norm)
  } else {
    Xb = Xnor
  }
  sparsity = sum(abs(htheta) > 0.001)
  sd.est = sqrt(sum((y - Xb %*% htheta)^2) / n)
  htheta = htheta * col.norm
  returnList = list("lasso.est" = htheta,
                     "sigma" = sd.est,
                     "sparsity" = sparsity)
  return(returnList)
}


### DLL: Original Version
### Function: Decorrelated Local Linear Estimator
# DLL = function(X, y, x.eval, f1.ind, h=NULL, data.swap=TRUE, quant.trans=TRUE, bwmethod="rot",nonpara=FALSE) {
#   n = nrow(X)
#   p = ncol(X)
#   n.eval = length(x.eval)
#   if (data.swap) {
#     a.ind = 1:round(n/2); b.ind = setdiff(1:n,a.ind)
#     n.a = length(a.ind); n.b = length(b.ind)
#   }
# 
#   # fit sparse additive model
#   if (data.swap) {
#     sam.a = cv.SAM(X[a.ind,],y[a.ind],quant.trans=quant.trans)
#     sam.b = cv.SAM(X[b.ind,],y[b.ind],quant.trans=quant.trans)
#     f.hat.a = predict.SAM(sam.b,X[a.ind,],quant.trans=quant.trans)
#     f.hat.b = predict.SAM(sam.a,X[b.ind,],quant.trans=quant.trans)
#     f.hat = rbind(f.hat.a,f.hat.b)
#   } else {
#     sam.model = cv.SAM(X,y,quant.trans=quant.trans)
#     f.hat = predict.SAM(sam.model,X,quant.trans=quant.trans)
#   }
# 
#   sigma1.sq = mean((y-apply(f.hat,1,sum))^2) # get sigma1.sq for variance estimation
# 
#   # calculate R by detracting nuisance function
#   R.hat = y - apply(f.hat[,-(f1.ind+1)], 1, sum) # f.ind plus 1 to skip the intercept
# 
# 
#   if (is.null(h)) {
#     if (bwmethod=="cv") {
#       h = suppressWarnings(regCVBwSelC(X[,f1.ind],R.hat,deg=1,kernel=SqK,interval=c(0.5,3)*n^(-0.2)))
#     }
#     if (bwmethod=="rot") {
#       h = suppressWarnings(thumbBw(X[,f1.ind],R.hat,deg=1,kernel=SqK))
#     }
#     if (bwmethod=="np") {
#       h = suppressWarnings(npregbw(R.hat~X[,f1.ind], ckertype="uniform",regtype = "ll"))$bw
#     }
#     if (bwmethod=="lp") {
#       h = lpbwselect(R.hat,X[,f1.ind],eval=x.eval,deriv=1,kernel="uni",bwselect="mse-rot")$bws[,"h"]
#     }
#   }
#   # try 
#   
#   # check for single bandwidth or a vector
#   if (length(h)==1) h = rep(h, n.eval)
# 
#   # calculate point estimator and standard error
#   est = est.se = Sn.hat = eff.samp = mean.delta = se.delta = rep(0, n.eval)
#   if (data.swap) {
#     # LASSO estimator to get delta.hat and mu.hat
#     if (!nonpara) {
#       gamma.a = Initialization.step(X[a.ind,-f1.ind], X[a.ind,f1.ind], lambda = "CV.min", intercept = TRUE)$lasso.est
#       gamma.b = Initialization.step(X[b.ind,-f1.ind], X[b.ind,f1.ind], lambda = "CV.min", intercept = TRUE)$lasso.est
#     } else {
#       sam.X1.a = cv.SAM(X[a.ind,-f1.ind],X[a.ind,f1.ind],quant.trans = TRUE)
#       sam.X1.b = cv.SAM(X[b.ind,-f1.ind],X[b.ind,f1.ind],quant.trans = TRUE)
#       X1.hat.a = predict(sam.X1.b$sam.final,newdata=X[a.ind,-f1.ind])$values
#       X1.hat.b = predict(sam.X1.a$sam.final,newdata=X[b.ind,-f1.ind])$values
#     }
# 
#     for (j in 1:n.eval) {
#       # calculate mu and delta
#       if (!nonpara) {
#         mu.a = x.eval[j] - cbind(1,matrix(X[a.ind,-f1.ind],n.a,p-1))%*%gamma.b
#         mu.b = x.eval[j] - cbind(1,matrix(X[b.ind,-f1.ind],n.b,p-1))%*%gamma.a
#         delta.a = matrix(X[a.ind,f1.ind],n.a,1) - cbind(1,matrix(X[a.ind,-f1.ind],n.a,p-1))%*%gamma.b
#         delta.b = matrix(X[b.ind,f1.ind],n.b,1) - cbind(1,matrix(X[b.ind,-f1.ind],n.b,p-1))%*%gamma.a
#         delta = c(delta.a,delta.b)
#         mean.delta[j] = mean(delta); se.delta[j] = sd(delta)
#       } else {
#         mu.a = x.eval[j] - X1.hat.a
#         mu.b = x.eval[j] - X1.hat.b
#         delta.a = matrix(X[a.ind,f1.ind],n.a,1) - X1.hat.a
#         delta.b = matrix(X[b.ind,f1.ind],n.b,1) - X1.hat.b
#         delta = c(delta.a,delta.b)
#         mean.delta[j] = mean(delta); se.delta[j] = sd(delta)
#       }
#       # calculate l
#       l.a = rep(0,n.a); l.b = rep(0,n.b)
#       for (i in 1:n.a) {
#         # weight
#         w.a = ifelse(abs(delta.a-mu.a[i])<=h[j],yes=1,no=0)
#         l.a[i] = sum((delta.a-mu.a[i])*w.a)/sum(w.a)
#         l.a[is.na(l.a)] = 0 # remove NAs
#       }
#       for (i in 1:n.b) {
#         # weight
#         w.b = ifelse(abs(delta.b-mu.b[i])<=h[j],yes=1,no=0)
#         l.b[i] = sum((delta.b-mu.b[i])*w.b)/sum(w.b)
#         l.b[is.na(l.b)] = 0 # remove NAs
#       }
#       l = c(l.a,l.b)
#       D.tilde = (X[,f1.ind]-x.eval[j]) - l
#       # uniform kernel
#       Kh = (1/h[j])*ifelse(abs(X[,f1.ind]-x.eval[j])/h[j]<=1, yes = 1/2, no = 0)
#       eff.samp[j] = sum(Kh)*2 # effective sample size in the bandwidth
#       D.hat = D.tilde - sum(D.tilde*Kh)/sum(Kh)
#       Sn.hat[j] = 1/n*sum(D.hat*(X[,f1.ind]-x.eval[j])*Kh)
#       est[j] = 1/(n*Sn.hat[j])*sum(D.hat*R.hat*Kh)
#       V.hat = sigma1.sq/(n^2*Sn.hat[j]^2)*sum(D.hat^2*Kh^2)
#       est.se[j] = sqrt(V.hat)
#     } # for each evaluation point
# 
# 
#   } else {
#     if (!nonpara) {
#       gamma = Initialization.step(X[,-f1.ind], X[,f1.ind], lambda = "CV.min", intercept = TRUE)$lasso.est
#     } else {
#       sam.X1 = cv.SAM(X[,-f1.ind],X[,f1.ind],quant.trans = T)
#       X1.hat = predict(sam.X1$sam.final,newdata=X[,-f1.ind])$values
#     }
#     for (j in 1:n.eval) {
#       # calculate mu and delta
#       if (!nonpara) {
#         mu = x.eval[j] - cbind(1,matrix(X[,-f1.ind],n,p-1))%*%gamma
#         delta= matrix(X[,f1.ind],n,1) - cbind(1,matrix(X[,-f1.ind],n,p-1))%*%gamma
#         mean.delta[j] = mean(delta); se.delta[j] = sd(delta)
#       } else {
#         mu = x.eval[j] - X1.hat
#         delta= matrix(X[,f1.ind],n,1) - X1.hat
#         mean.delta[j] = mean(delta); se.delta[j] = sd(delta)
#       }
#       # calculate l
#       l = rep(0,n)
#       for (i in 1:n) {
#         # weight
#         w = ifelse(abs(delta-mu[i])<=h[j],yes=1,no=0)
#         l[i] = sum((delta-mu[i])*w)/sum(w)
#         l[is.na(l)] = 0 # remove NAs
#       }
#       D.tilde = (X[,f1.ind]-x.eval[j]) - l
#       # uniform kernel
#       Kh = (1/h[j])*ifelse(abs(X[,f1.ind]-x.eval[j])/h[j]<=1, yes = 1/2, no = 0)
#       eff.samp[j] = sum(Kh)*2 # effective sample size in the bandwidth
#       D.hat = D.tilde - sum(D.tilde*Kh)/sum(Kh)
#       Sn.hat[j] = 1/n*sum(D.hat*(X[,f1.ind]-x.eval[j])*Kh)
#       est[j] = 1/(n*Sn.hat[j])*sum(D.hat*R.hat*Kh)
#       V.hat = sigma1.sq/(n^2*Sn.hat[j]^2)*sum(D.hat^2*Kh^2)
#       est.se[j] = sqrt(V.hat)
#     } # for each evaluation point
#   }
# 
#   if (data.swap) {
#     returnList = list(est = est,
#                        est.se = est.se,
#                        x.eval = x.eval,
#                        h = h,
#                        Sn.hat = Sn.hat,
#                        eff.samp = eff.samp,
#                        R.hat = R.hat,
#                        f.hat = f.hat,
#                        sam.a = sam.a, sam.b = sam.b,
#                        sigma1.sq = sigma1.sq,
#                        mean.delta = mean.delta,
#                        se.delta = se.delta,
#                        gamma.a = gamma.a, gamma.b = gamma.b,
#                        l.a = l.a, l.b = l.b
#                        )
#   } else {
#     returnList = list(est = est,
#                        est.se = est.se,
#                        x.eval = x.eval,
#                        h = h,
#                        Sn.hat = Sn.hat,
#                        eff.samp = eff.samp,
#                        R.hat = R.hat,
#                        f.hat = f.hat,
#                        sam.model = sam.model,
#                        sigma1.sq = sigma1.sq,
#                        mean.delta = mean.delta,
#                        se.delta = se.delta,
#                        gamma = gamma,
#                        l = l
#                        )
#   }
#   returnList
# }


### DLL
### Function: Decorrelated Local Linear Estimator
### Extended: foi.ind can be a vector now
DLL = function(X, y, x.eval, foi.ind, h=NULL, data.swap=TRUE, quant.trans=TRUE, bwmethod="rot",nonpara=FALSE) {
  n = nrow(X)
  p = ncol(X)
  n.eval = length(x.eval)
  if (data.swap) {
    a.ind = 1:round(n/2); b.ind = setdiff(1:n,a.ind)
    n.a = length(a.ind); n.b = length(b.ind)
  }
  
  # fit sparse additive model
  if (data.swap) {
    sam.a = cv.SAM(X[a.ind,],y[a.ind],quant.trans=quant.trans)
    sam.b = cv.SAM(X[b.ind,],y[b.ind],quant.trans=quant.trans)
    f.hat.a = predict.SAM(sam.b,X[a.ind,],quant.trans=quant.trans)
    f.hat.b = predict.SAM(sam.a,X[b.ind,],quant.trans=quant.trans)
    f.hat = rbind(f.hat.a,f.hat.b)
  } else {
    sam.model = cv.SAM(X,y,quant.trans=quant.trans)
    f.hat = predict.SAM(sam.model,X,quant.trans=quant.trans)
  }
  
  sigma1.sq = mean((y-apply(f.hat,1,sum))^2) # get sigma1.sq for variance estimation
  
  # calculate point estimator and standard error
  est = est.se = matrix(NA,length(x.eval),length(foi.ind))
  colnames(est) = colnames(est.se) = paste("f",as.character(foi.ind),sep="")
  rownames(est) = rownames(est.se) = as.character(x.eval)
  bw.save = rep(NA,length(foi.ind))
  
  for (f1.ind in 1:length(foi.ind)) {
    # calculate R by detracting nuisance function
    R.hat = y - apply(f.hat[,-(foi.ind[f1.ind]+1)], 1, sum) # foi.ind[f1.ind] plus 1 to skip the intercept
    
    
    if (bwmethod=="cv") {
      h = suppressWarnings(regCVBwSelC(X[,foi.ind[f1.ind]],R.hat,deg=1,kernel=SqK,interval=c(0.5,3)*n^(-0.2)))
    }
    if (bwmethod=="rot") {
      h = suppressWarnings(thumbBw(X[,foi.ind[f1.ind]],R.hat,deg=1,kernel=SqK))
    }
    if (bwmethod=="np") {
      h = suppressWarnings(npregbw(R.hat~X[,foi.ind[f1.ind]], ckertype="uniform",regtype = "ll"))$bw
    }
    if (bwmethod=="lp") {
      h = lpbwselect(R.hat,X[,foi.ind[f1.ind]],eval=x.eval,deriv=1,kernel="uni",bwselect="mse-rot")$bws[,"h"]
    }
    bw.save[f1.ind] = h
    # check for single bandwidth or a vector
    if (length(h)==1) h = rep(h, n.eval)
    
    if (data.swap) {
      # LASSO estimator to get delta.hat and mu.hat
      if (!nonpara) {
        gamma.a = Initialization.step(X[a.ind,-foi.ind[f1.ind]], X[a.ind,foi.ind[f1.ind]], lambda = "CV.min", intercept = TRUE)$lasso.est
        gamma.b = Initialization.step(X[b.ind,-foi.ind[f1.ind]], X[b.ind,foi.ind[f1.ind]], lambda = "CV.min", intercept = TRUE)$lasso.est
      } else {
        sam.X1.a = cv.SAM(X[a.ind,-foi.ind[f1.ind]],X[a.ind,foi.ind[f1.ind]],quant.trans = TRUE)
        sam.X1.b = cv.SAM(X[b.ind,-foi.ind[f1.ind]],X[b.ind,foi.ind[f1.ind]],quant.trans = TRUE)
        X1.hat.a = predict(sam.X1.b$sam.final,newdata=X[a.ind,-foi.ind[f1.ind]])$values
        X1.hat.b = predict(sam.X1.a$sam.final,newdata=X[b.ind,-foi.ind[f1.ind]])$values
      }
      
      for (j in 1:n.eval) {
        # calculate mu and delta
        if (!nonpara) {
          mu.a = x.eval[j] - cbind(1,matrix(X[a.ind,-foi.ind[f1.ind]],n.a,p-1))%*%gamma.b
          mu.b = x.eval[j] - cbind(1,matrix(X[b.ind,-foi.ind[f1.ind]],n.b,p-1))%*%gamma.a
          delta.a = matrix(X[a.ind,foi.ind[f1.ind]],n.a,1) - cbind(1,matrix(X[a.ind,-foi.ind[f1.ind]],n.a,p-1))%*%gamma.b
          delta.b = matrix(X[b.ind,foi.ind[f1.ind]],n.b,1) - cbind(1,matrix(X[b.ind,-foi.ind[f1.ind]],n.b,p-1))%*%gamma.a
          delta = c(delta.a,delta.b)
        } else {
          mu.a = x.eval[j] - X1.hat.a
          mu.b = x.eval[j] - X1.hat.b
          delta.a = matrix(X[a.ind,foi.ind[f1.ind]],n.a,1) - X1.hat.a
          delta.b = matrix(X[b.ind,foi.ind[f1.ind]],n.b,1) - X1.hat.b
          delta = c(delta.a,delta.b)
        }
        # calculate l
        l.a = rep(0,n.a); l.b = rep(0,n.b)
        for (i in 1:n.a) {
          # weight
          w.a = ifelse(abs(delta.a-mu.a[i])<=h[j],yes=1,no=0)
          l.a[i] = sum((delta.a-mu.a[i])*w.a)/sum(w.a)
          l.a[is.na(l.a)] = 0 # remove NAs
        }
        for (i in 1:n.b) {
          # weight
          w.b = ifelse(abs(delta.b-mu.b[i])<=h[j],yes=1,no=0)
          l.b[i] = sum((delta.b-mu.b[i])*w.b)/sum(w.b)
          l.b[is.na(l.b)] = 0 # remove NAs
        }
        l = c(l.a,l.b)
        D.tilde = (X[,foi.ind[f1.ind]]-x.eval[j]) - l
        # uniform kernel
        Kh = (1/h[j])*ifelse(abs(X[,foi.ind[f1.ind]]-x.eval[j])/h[j]<=1, yes = 1/2, no = 0)
        D.hat = D.tilde - sum(D.tilde*Kh)/sum(Kh)
        Sn.hat = 1/n*sum(D.hat*(X[,foi.ind[f1.ind]]-x.eval[j])*Kh)
        est[j,f1.ind] = 1/(n*Sn.hat)*sum(D.hat*R.hat*Kh)
        V.hat = sigma1.sq/(n^2*Sn.hat^2)*sum(D.hat^2*Kh^2)
        est.se[j,f1.ind] = sqrt(V.hat)
      } # for each evaluation point
      
      
    } else {
      if (!nonpara) {
        gamma = Initialization.step(X[,-foi.ind[f1.ind]], X[,foi.ind[f1.ind]], lambda = "CV.min", intercept = TRUE)$lasso.est
      } else {
        sam.X1 = cv.SAM(X[,-foi.ind[f1.ind]],X[,foi.ind[f1.ind]],quant.trans = T)
        X1.hat = predict(sam.X1$sam.final,newdata=X[,-foi.ind[f1.ind]])$values
      }
      for (j in 1:n.eval) {
        # calculate mu and delta
        if (!nonpara) {
          mu = x.eval[j] - cbind(1,matrix(X[,-foi.ind[f1.ind]],n,p-1))%*%gamma
          delta= matrix(X[,foi.ind[f1.ind]],n,1) - cbind(1,matrix(X[,-foi.ind[f1.ind]],n,p-1))%*%gamma
        } else {
          mu = x.eval[j] - X1.hat
          delta = matrix(X[,foi.ind[f1.ind]],n,1) - X1.hat
        }
        # calculate l
        l = rep(0,n)
        for (i in 1:n) {
          # weight
          w = ifelse(abs(delta-mu[i])<=h[j],yes=1,no=0)
          l[i] = sum((delta-mu[i])*w)/sum(w)
          l[is.na(l)] = 0 # remove NAs
        }
        D.tilde = (X[,foi.ind[f1.ind]]-x.eval[j]) - l
        # uniform kernel
        Kh = (1/h[j])*ifelse(abs(X[,foi.ind[f1.ind]]-x.eval[j])/h[j]<=1, yes = 1/2, no = 0)
        D.hat = D.tilde - sum(D.tilde*Kh)/sum(Kh)
        Sn.hat = 1/n*sum(D.hat*(X[,foi.ind[f1.ind]]-x.eval[j])*Kh)
        est[j,f1.ind] = 1/(n*Sn.hat)*sum(D.hat*R.hat*Kh)
        V.hat = sigma1.sq/(n^2*Sn.hat^2)*sum(D.hat^2*Kh^2)
        est.se[j,f1.ind] = sqrt(V.hat)
      } # for each evaluation point
    }
  }
 
  
  returnList = list(est = est,
                    est.se = est.se,
                    x.eval = x.eval,
                    bw.save = bw.save,
                    sigma1.sq = sigma1.sq
  )
  
  returnList
}

