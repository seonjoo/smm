#' Conduct K-fold cross validation for sparse moderated mediation with group lasso
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param K (default=5) number of cross-validation folds
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:30)/100)) tuning parameter for L1 penalization
#' @param grpgroup (default=c(rep(1,3),rep( 1:V +1,5)))
#' @param penalty.factor (default=c(0,rep(sqrt(2),V))) give different weight of penalization for the 2V mediation paths.
#' @param multicore (default=1) number of multicore
#' @param seednum (default=10000) seed number for cross validation
#' @return cv.lambda: optimal lambda
#' @return cv.mse: minimum MSE value
#' @return mse: Array of MSE, length(alpha) x length(lambda) x length (tau)
#' @return lambda: vector of lambda
#' @return tau: vector of tau used
#' @return alpha: vector of alpha used
#' @return z: cross-valication results
#' @examples
#' N=200
#' V=50
#' set.seed(1234)
#' a = rbinom(V,1,0.1)*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  as.vector(X + M %*% b + rnorm(N))
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords hdlfpca gglasso
#' @import parallel
#' @import MASS
#' @import gglasso
#' @export


cv.sparse.modmediation.grplasso= function(X,M,Y,Z,tol=10^(-5),K=5,max.iter=100,
                                       lambda= log(1+(1:15)/40),
                                       grpgroup=c(rep(1,3), rep(1:ncol(M)+1,5)),
                                       penalty.factor=c(0,rep(1,ncol(M))),
                                       verbose=FALSE,
                                       multicore=1,seednum=1000000){
  ## Center all values
  N = nrow(M)
  V = ncol(M)
  #Y.mean=mean(Y)
  #X.mean=mean(X)
  #M.mean=apply(M,2,mean)
  #Y.sd=sqrt(var(Y))
  #X.sd=sqrt(var(X))
  #M.sd=sqrt(apply(M,2,var))

  Y = scale(Y,center=TRUE,scale=TRUE)
  X = matrix(scale(X,center=TRUE,scale=TRUE),N,1)
  M = scale(M, center=TRUE,scale=TRUE)
  Z = scale(Z, center=TRUE,scale=TRUE)

  ###K-fold cross-validation
  set.seed(seednum)
  cvid = (rep(1:K, each=ceiling(N/K))[1:N])[sort.int(rnorm(N),index.return=TRUE)$ix]
##  sparse.mediation.grplasso.fold(1, Y,X,M,cvid,lambdas[10], max.iter, tol)
  if(multicore>1){
    options(cores = multicore)
    z<-mclapply(1:K, function(fold){
      sparse.modmediation.grplasso.fold(fold,X,M,Y,Z,cvid,tol=tol,max.iter=max.iter,
                                              lambda = lambda,
                                              grpgroup=grpgroup,
                                              penalty.factor=penalty.factor,
                                              threshold=threshold,
                                              verbose=verbose)
#    sparse.mediation.grplasso.fold(fold, Y,X,M,cvid,lambda, max.iter, tol)
      }, mc.cores=multicore)
  }else{
    z<-lapply(1:K, function(fold){sparse.modmediation.grplasso.fold(fold,X,M,Y,Z,cvid,tol=tol,max.iter=max.iter,
                                                                          lambda = lambda,
                                                                          grpgroup=grpgroup,
                                                                          penalty.factor=penalty.factor,
                                                                          threshold=threshold,
                                                                          verbose=verbose)})
  }

  mseest=apply(do.call(cbind,lapply(z,function(x)x$mse$mse)),1,sum)
  minloc=which.min(mseest)
  min.lambda=lambda[minloc]


  return(list(cv.lambda=min.lambda, cv.mse=mseest[minloc],
              mse=mseest, lambda=lambda,z=z))

}

sparse.modmediation.grplasso.fold<-function(fold, X,M,Y,Z,cvid,
                                         tol=10^(-10),max.iter=100,
                                         lambda=c(0.2,0.4),
                                         grpgroup=c(rep(1,3), rep(1:V+1,5)),
                                         penalty.factor=c(0,rep(1,V)),
                                         threshold=0.00001,
                                         verbose=FALSE){
  test.indx=which(cvid==fold)
  train.indx=which(cvid!=fold)
  fit = sparse.modmediation.grplasso(X[train.indx,],M[train.indx,],Y[train.indx,],Z[train.indx,],
                                          tol=tol, max.iter=max.iter,
                                    lambda =lambda,
                                    grpgroup=grpgroup,
                                    penalty.factor=penalty.factor,
                                    threshold=threshold,
                                    verbose=verbose)
  fit$mse = cv.sparse.modmediation.msecomputing(obj=fit, Y.test=Y[test.indx,], X.test=X[test.indx,],
                                                M.test=M[test.indx,], Z.test=Z[test.indx,])
  return(fit)
}

cv.sparse.modmediation.msecomputing<-function(obj, Y.test, X.test, M.test, Z.test){
  #	print(obj)
  yhat = cbind(X.test, X.test*Z.test, Z.test) %*% obj$c + M.test %*% obj$hatb1 + (as.vector(Z.test)*M.test) %*% obj$hatb2
#  mhat = X.test %*% obj$hata1  + (Z.test*X.test) %*% obj$hata1
  mse.m = rep(0,length(obj$lambda))
  for (j in 1:length(obj$lambda)){
    mhat=X.test %*% t( obj$hata1[,j]) + (Z.test*X.test) %*% t(obj$hata2[,j])
    mse.m[j]=sum((M.test - mhat)^2)
  }
  mse=(apply(Y.test-yhat, 2,function(x){sum(x^2)})) + mse.m
  return(list(mse=mse,lambda=obj$lambda))
}
