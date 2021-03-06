#' Conduct sparse mediation with group LASSO
#'
#' Fit a mediation model via penalized maximum likelihood and structural equation model.
#' The regularization path is computed using group lasso at a grid of
#' values for the regularization parameter lambda. Currently, mediation analysis is developed
#' based on gaussian assumption.
#'
#' Multiple Mediaton Model:
#' (1) M = Xa + e1
#' (2) Y = Xc' + Mb + e2
#' And in the optimization, we do not regularize c', due to the assumption of partial mediation.
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda1 (default=exp(-5:0)) tuning parameter for c',a, b coefficients
#' @param lambda2 (default=exp(seq(0,0.5*log(ncol(M)),length=3))) tuning parameter for the Omega=Sigma_2^{-1}
#' @param alpha (default=0) alpha=0, group lasso penalization will run
#' @param group.penalty.factor (V+1)-dimensional group penalty factor vector. If a user does not want to penalize mediator, specify 0 otherwise 1. The first element is the direct effect followed by V-mediators. The default value is c(0,rep(1,V)).
#' @param penalty.factor (1+2*V)-dimensional sparsity penalty factor vector.
#' @param threshold (default=10^(-8))
#' @param non.zeros.stop (default=ncol(M)) when to stop searching regularization path
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
#' a = c(rep(1,3),rep(0,V-3));b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' fit=sparse.mediation.grplasso(X,M,Y,verbose=FALSE,
#' lambda1 = exp(seq(-5,-1, length=10)),
#' alpha=c(0.95,0.75,0.5),lambda2=c(1,2),Omega.out = FALSE)
#' fit
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation glmnet
#' @import parallel
#' @import MASS
#' @import gglasso
#' @importFrom stats var predict
#' @export
sparse.mediation.grplasso = function(X,M,Y,
                                     tol=10^(-10),
                                     max.iter=100,
                                     lambda1 = 0.01,#exp(-5:0),
                                     lambda2 = 1,#exp(seq(0,0.5*log(ncol(M)),length=3)),
                                     alpha=0,
                                     group.penalty.factor=c(1, rep(1, ncol(M))),
                                     penalty.factor=c(1,rep(1,ncol(M))),
                                     Omega.out=FALSE,
                                     verbose=FALSE,
                                     threshold=0,
                                     non.zeros.stop=ncol(M)){

    re=list()
    V = ncol(M)
    N=nrow(M)

    # if( V< min(200, floor(N*0.75))){
    #   re=sparse.mediation.grplasso.smallp(X,M,Y,tol=tol,max.iter=max.iter,
    #                                       lambda = lambda1,verbose=verbose,
    #                                       penalty.factor=penalty.factor,
    #                                       grpgroup=grpgroup)
    #  }else{
       # re=sparse.mediation.grplasso.largep_omega(X,M,Y,tol=tol,max.iter=max.iter,
       #                                           lambda1 = lambda1,
       #                                           lambda2 = lambda2,
       #                                           alpha=alpha,
       #                                           verbose=verbose,
       #                                           Omega.out=Omega.out,
       #                                           penalty.factor=penalty.factor,
       #                                           grpgroup=grpgroup)

       re=sparse.mediation.sgrplasso.largep_omega(X,M,Y,tol=tol,max.iter=max.iter,
                                                 lambda1 = lambda1,
                                                 lambda2 = lambda2,
                                                 alpha=alpha,
                                                 verbose=verbose,
                                                 Omega.out=Omega.out,
                                                 penalty.factor=penalty.factor,
                                                 group.penalty.factor=group.penalty.factor,
                                                 threshold=threshold,
                                                 non.zeros.stop=non.zeros.stop
                                                 )


           # }
    return(re)
}
