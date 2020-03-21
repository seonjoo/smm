#' Conduct K-fold cross validation for sparse mediation with group lasso with multiple tuning parameters
#'
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param K (default=5) number of cross-validation folds
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda1 (default=exp(-5:0)) tuning parameter for L1 penalization
#' @param lambda2 (default=exp(-5:0)) tuning parameter for L1 penalization
#' @param penalty.factor (default=c(0,rep(1,2*V))) give different weight of penalization for the 2V mediation paths.
#' @param multicore (default=1) number of multicore
#' @param seednum (default=10000) seed number for cross validation
#' @return cv.lambda: optimal lambda
#' @return cv.tau: optimal tau
#' @return cv.alpha: optimal tau
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
#' a = c(rep(1,3),rep(0,V-3))*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  as.vector(X + M %*% b + rnorm(N))
#' system.time(cvfit<-cv.sparse.mediation.grplasso(X, M, Y, K = 4,multicore = 4, seednum = 20200317))
#' cvfit$cv.lambda
#' fit<-sparse.mediation.grplasso(X,M,Y,lambda1 = cvfit$cv.lambda1,lambda2 = cvfit$cv.lambda2)
#' nonzerogroups = 1-as.numeric((fit$hata!=0)+(fit$hatb!=0) ==0)
#' refit<-sparse.mediation.grplasso(X,M[,nonzerogroups==1],Y,lambda = 0,lambda2 = cvfit$cv.lambda2)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords gglasso
#' @import parallel
#' @import MASS
#' @import gglasso
#' @import Matrix
#' @export


cv.sparse.mediation.grplasso= function(X,M,Y,tol=10^(-5),K=5,max.iter=100,
                                       lambda1= exp(-5:0),
                                       lambda2= exp(seq(0,0.5*log(ncol(M)),length=3)),
                                       grpgroup=c(1, rep(1:(ncol(M))+1,2)),
                                       penalty.factor=c(0,rep(1,ncol(M))),
                                       verbose=FALSE,
                                       multicore=1,seednum=100000){
  ## Center all values
  N = nrow(M)
  V = ncol(M)
  Y.mean=mean(Y)
  X.mean=mean(X)
  M.mean=apply(M,2,mean)
  Y.sd=sqrt(var(Y))
  X.sd=sqrt(var(X))
  M.sd=sqrt(apply(M,2,var))

  Y = scale(Y,center=TRUE,scale=TRUE)
  X = matrix(scale(X,center=TRUE,scale=TRUE),N,1)
  M = scale(M, center=TRUE,scale=TRUE)

  ###K-fold cross-validation
  set.seed(seednum)
  cvid = (rep(1:K, each=ceiling(N/K))[1:N])[sort.int(rnorm(N),index.return=TRUE)$ix]
##  sparse.mediation.grplasso.fold(1, Y,X,M,cvid,lambdas[10], max.iter, tol)
  if(multicore>1){
#    options(cores = multicore)
    z<-mclapply(1:K, function(fold){
      sparse.mediation.grplasso.fold(fold, Y,X,M,cvid,lambda1=lambda1, lambda2=lambda2, max.iter=max.iter, tol=tol)}, mc.cores=multicore)
  }else{
    z<-lapply(1:K, function(fold){sparse.mediation.grplasso.fold(fold, Y,X,M,cvid,lambda1, lambda2,max.iter, tol)})
  }

  mseest=apply(do.call(cbind,lapply(z,function(x)x$mse$mse)),1,sum)
  minloc=which.min(mseest)
  min.lambda1=z[[1]]$lambda1[minloc]
  min.lambda2=z[[1]]$lambda2[minloc]


  return(list(cv.lambda1=min.lambda1, cv.lambda2=min.lambda2,cv.mse=mseest[minloc],
              mse=mseest, lambda1=z[[1]]$lambda1,z=z, lambda2=z[[1]]$lambda2))

}



