#' Compare x to 1
#' @param obj
#' @param Y.test
#' @param X.test
#' @param M.test
#' @keywords internal
#' @export
#' @NoRd
cv.sparse.mediation.msecomputing<-function(obj, Y.test, X.test, M.test){
  V = (nrow(obj$beta))
  a.train=obj$hata
  b.train=obj$hatb
  c.train=matrix(obj$c,nrow=1)
  #	print(obj)
  yhat = X.test %*% c.train + M.test %*% b.train
  mse.m = rep(0,length(obj$lambda1))
  for (j in 1:length(obj$lambda1)){
    mhat=X.test %*% t(a.train[,j])
    mse.m[j]=sum((M.test - mhat)^2)
  }
  mse=(apply(Y.test-yhat, 2,function(x){sum(x^2)})) + mse.m
  return(list(mse=mse,lambda1=obj$lambda1,lambda2=obj$lambda2))
}
