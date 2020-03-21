#' Compute BIC scores for model selection
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:50)/125)) tuning parameter for L1 penalization
#' @param grpgroup (default=c(1,rep( 1:V +1,2)))
#' @param penalty.factor (default=c(0,rep(sqrt(2),V))) give different weight of penalization for the 2V mediation paths.
#' @return c directeffect
#' @return hatb Path b (M->Y given X) estimates
#' @return hata Path a (X->M) estimates
#' @return medest Mediation estimates (a*b)
#' @return alpha
#' @return lambda
#' @return nump Number of selected mediation paths
#' @examples
#' N=100
#' V=50
#' set.seed(1234)
#' a = rbinom(V,1,0.1)*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation glmnet
#' @import parallel
#' @import MASS
#' @import gglasso
#' @export
#'
#'
computebic<-function(X,Y,M,a,b,c,tol=10^(-10),
                     max.iter=100,grpgroup=c(1, rep(1:V+1,2)),
                     threshold=0){
  N=nrow(M);V=ncol(M)
  zerolist=c(0,as.numeric(b*a==0))
  refit<-NULL
  try(refit<-sparse.mediation.grplasso(X,M,Y,tol=tol,max.iter=max.iter,
                                             lambda = 10000,
                                             grpgroup=grpgroup,
                                             penalty.factor=zerolist,
                                             threshold=threshold,
                                             verbose=FALSE))

      if(is.null(refit)==FALSE){
    tmp = M - matrix(X, N, 1) %*% matrix(refit$hata, 1,V)
    Sigma2 = t(tmp) %*% tmp/N
    bic=N*log(sum(Y - cbind(X,M) %*% c(refit$c,refit$hatb))^2/N) +N*log(det(Sigma2)) + log(N)*(sum(1-zerolist))
  }
if(is.null(refit)==TRUE){
  refit=list(medest=rep(NA,V))
  bic=Inf
}
return(list(bic=bic, refit=refit))# + V*(1+V)/2+1))
}
