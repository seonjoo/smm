#' sparse.mediation.grplasso.fold
#' @param fold
#' @keywords internal
#' @export
#' @NoRd

sparse.mediation.grplasso.fold<-function(fold, Y,X,M,cvid,
                                         lambda1=0.01,
                                         lambda2=1,
                                         tol=10^(-10),max.iter=100,
                                         grpgroup=c(1, rep(1:(ncol(M))+1,2)),
                                         penalty.factor=c(0,rep(1,ncol(M))),
                                         threshold=0.00001,
                                         verbose=FALSE){
  test.indx=which(cvid==fold)
  train.indx=which(cvid!=fold)
  fit = sparse.mediation.grplasso(Y=Y[train.indx], X=X[train.indx,], M=M[train.indx,],
                                  lambda1=lambda1,
                                  lambda2=lambda2,
                                  tol=tol,
                                  max.iter=max.iter,
                                  grpgroup=grpgroup,
                                  penalty.factor=penalty.factor,
                                  verbose=verbose )

  fit$mse = cv.sparse.mediation.msecomputing(obj=fit, Y.test=Y[test.indx], X.test=X[test.indx,], M.test=M[test.indx,])
  return(fit)
}
