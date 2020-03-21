#' Conduct sparse mediation with group LASSO for smaller p (< min(n,200))
#'
#' Fit a mediation model via penalized maximum likelihood and structural equation model.
#' The regularization path is computed for the lasso or elasticnet penalty at a grid of
#' values for the regularization parameter lambda. Currently, mediation analysis is developed based on gaussian assumption.
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
#' @param lambda (default=exp(-5:0)) tuning parameter for L1 penalization
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
#' a = c(rep(1,3),rep(0,V-3))*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' fit=sparse.mediation.grplasso.smallp(X,M,Y,verbose=FALSE, lambda = exp(-5:0))
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation glmnet
#' @import parallel
#' @import MASS
#' @import gglasso
#' @importFrom stats var predict
#' @export
sparse.mediation.grplasso.smallp = function(X,M,Y,
                                     tol=10^(-10),
                                     max.iter=100,
                                     lambda = exp(-5:0),
                                     grpgroup=c(1, rep(1:(ncol(M))+1,2)),
                                     penalty.factor=c(0,rep(1,ncol(M))),
                                     verbose=FALSE){


  ## Center all values, and also make their scales to be 1. In this context, all coefficients will be dexribed in terms of correlation or partial correlations.
  N = nrow(M)
  V = ncol(M)
  Y.mean=mean(Y)
  X.mean=mean(X)
  M.mean=apply(M,2,mean)
  Y.sd=as.vector(sqrt(var(Y)))
  X.sd=as.vector(sqrt(var(X)))
  M.sd=sqrt(apply(M,2,var))
  Y = scale(Y,center=TRUE,scale=TRUE)
  X = matrix(scale(X,center=TRUE,scale=TRUE),N,1)
  M = scale(M, center=TRUE,scale=TRUE)

  ## Penalty Factor
  if (ncol(X)>1){stop("X has more than 1 colum. Stop.")}
  ## Initialization###
  ## OLS Estimation ###
  U = cbind(X,M)
  tUU = t(U)%*%U
  tUU.sqmat=sqrtmat.comp(tUU)
  invtUU = ginv(tUU)
  #invtMM = ginv(t(M)%*%M)
  tXX = t(X)%*%X
  tUY = t(U)%*%Y
  tMX = t(M)%*%X

  ## Interative Update
  lam=lambda
  betaest =  matrix(0,1+2*V,length(lam))

  for (j in 1:length(lam)){
    if (verbose==TRUE){print(paste("Lambda",lam[j]))}
    gamma_new = invtUU %*% tUY
    alpha_new = t(ginv(t(X)%*%X)%*%t(X)%*%M)

    iter=0
    err=1000

    while( err>tol & iter<max.iter){

      alpha_old=alpha_new
      gamma_old = gamma_new
      beta_old = c(gamma_old,alpha_old)

      sigma1 = mean((Y - U %*% gamma_old)^2)
      tmp = M - matrix(X,N,1) %*% matrix(alpha_old,1,V)
      Sigma2 = t(tmp)%*%tmp/N
      Sigma2.sqmat=sqrtmat.comp(Sigma2)
      Sigma2.sqrt.inv=ginv(Sigma2.sqmat)
      Sigma2.inv=ginv(Sigma2)

      A = matrix(0,1+2*V,1+2*V)
      A[1:(1+V),1:(1+V)]=1/sigma1 * tUU
      A[(1+V)+ 1:V,(1+V)+ 1:V]=  as.numeric(tXX) * Sigma2.inv

      sqmatA = A;sqmatA[1:(1+V),1:(1+V)]=1/sqrt(sigma1) * tUU.sqmat
      sqmatA[(1+V)+ 1:V,(1+V)+ 1:V]=  sqrt(as.numeric(tXX)) * Sigma2.sqrt.inv
      C = ginv(sqmatA) %*% rbind(tUY/sigma1, Sigma2.inv%*%tMX)

      if(is.null(penalty.factor)==TRUE){
        #fit = glmnet(sqmatA, C,lambda=lambda[j],alpha=alpha)
        fit=gglasso(x=scale(sqmatA)[,order(grpgroup)],
                    y=scale(C),
                    lambda=lambda[j],
                    group=grpgroup[order(grpgroup)])
      }else{
        #fit = glmnet(sqmatA, C,lambda=lambda[j],penalty.factor=penalty.factor,alpha=alpha)
        fit=gglasso(x=scale(sqmatA)[,order(grpgroup)],
                    y=scale(C),
                    lambda=lambda[j],
                    group=grpgroup[order(grpgroup)],
                    pf = penalty.factor)
      }

      beta_new = as.vector(coef(fit))[-1]
      ## use thresholds as well: since all variables are standardized, coefficients less than 0.001 does not have any meaning.
      #if (threshold>0){
      #  beta_new[abs(beta_new)<threshold]<-0
      #}
      #beta_new[(1:V) +1]*beta_new[(1:V) +V+1]
      gamma_new = beta_new[c(1, (1:V)*2)]#beta_new[1:(V+1)]
      alpha_new = beta_new[c(1:V)*2+1]#beta_new[(1:V)+ V+1]
      err = sum((beta_old[-1]-c(gamma_new[-1],alpha_new))^2)
      iter=iter+1
      if (verbose==TRUE){print(c(iter, err,sigma1))}
    }
    betaest[,j]=beta_new[c(1, (1:V)*2, (1:V)*2+1)]
  }
  cest =betaest[1,]
  medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]
  nump=apply(betaest,2,function(x){sum(abs(x)>0)})


  return(list(
    c = cest,
    hatb=betaest[(1:V)+1,]*Y.sd/M.sd,
    hata=betaest[(1:V)+V+1,]*M.sd/X.sd,
    medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]*Y.sd/X.sd,
    lambda = lambda,
    nump=nump
  ))
}
