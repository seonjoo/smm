#' Conduct sparse mediation for large p ( p > n)
#'
#' Sparse mediation with sparse group lasso (for mediation paths)
#' and sparse precision matrix estimation using fast computation of inverse matrix
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
#' @param lambda2 Tuning parameter for Covariance matrix L1 penalization
#' @param lambda1 (default=seq(0.02,0.4,length=5)) tuning parameter for regression coefficient L1 penalization
#' @param alpha (defult=1) tuning parameter for L2 penalization
#' @param verbose (default=FALSE) print progress
#' @param Omega.out (defult=FALSE) output Omega estimates (beta version WIP.)
#' @param threshold (default=10^(-5))
#' @param non.zeros.stop (default=ncol(M)) When to stop the regularization path.
#' @param group.penalty.factor (V+1)-dimensional group penalty factor vector. If a user does not want to penalize mediator, specify 0 otherwise 1. The first element is the direct effect followed by V-mediators. The default value is c(0,rep(1,V)).
#' @param penalty.factor (1+2*V)-dimensional sparsity penalty factor vector.
#' @return c directeffect
#' @return hatb Path b (M->Y given X) estimates
#' @return hata Path a (X->M) estimates
#' @return medest Mediation estimates (a*b)
#' @return alpha
#' @return lambda1 Tuning parameters for regression coefficients
#' @return lambda2 Tuning parameters for inversed covariance matrix (Omega)
#' @return nump Number of selected mediation paths
#' @return Omega Estimated covariance matrix of the mediator
#' @examples
#' library(sparsemediation)
#' N=100
#' V=50
#' set.seed(1234)
#' covmat=matrix(0,V+2,V+2);
#' covmat[1,2]=0.5;covmat[1, (1:3)+2]=rep(0.5,3);covmat[2, (1:3)+2]=rep(0.5,3);
#' covmat=covmat+t(covmat);diag(covmat)<-1
#' sqrtmat = sqrtmat.comp(covmat)
#' tmpmat = matrix(rnorm(N*(V+2)),N,V+2) %*% sqrtmat
#'
#' X=tmpmat[,1]
#' Y=tmpmat[,2]
#' M=tmpmat[,-c(1:2)]
#' #fit=sparse.mediation.sgrlasso.largep_omega(X,M,Y)
#'
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation L1penalization
#' @import parallel
#' @import MASS
#' @import glmnet
#' @import QUIC
#' @import Matrix
#' @importFrom stats var predict
#' @useDynLib smm
#' @importFrom Rcpp sourceCpp
#' @export
sparse.mediation.sgrplasso.largep_omega = function(X,M,Y,#Cov=NULL,
                                         tol=10^(-10),
                                         max.iter=10,
                                         lambda1 =exp(-5:1),
                                         lambda2 = exp(seq(0,0.8*log(ncol(M)),length=4)),
                                         alpha=0.8,
                                         group.penalty.factor=c(1,rep(1, ncol(M))),
                                         penalty.factor=c(1,rep(1,ncol(M)*2)),
                                         verbose=FALSE,
                                         Omega.out=FALSE,
                                         threshold=10^(-5),
                                         non.zeros.stop=ncol(M)/2){


  ## Center all values, and also make their scales to be 1. In this context, all coefficients will be dexribed in terms of correlation or partial correlations.
  N = nrow(M)
  V = ncol(M)
  #Y.mean=mean(Y)
  #X.mean=mean(X)
  #M.mean=apply(M,2,mean)
  Y.sd=as.vector(sqrt(var(Y)))
  X.sd=as.vector(sqrt(var(X)))
  M.sd=sqrt(apply(M,2,var))
  Y = scale(Y,center=FALSE,scale=TRUE)
  X = matrix(scale(X,center=TRUE,scale=TRUE),N,1)
  M = scale(M, center=FALSE,scale=TRUE)

#  if(is.null(Cov)==FALSE){Cov = scale(Cov, center=FALSE,scale=TRUE)}

  ## Penalty Factor
  if (ncol(X)>1){stop("X has more than 1 colum. Stop.")}
  ## Initialization###
  ## OLS Estimation ###
  U = cbind(X,M)

  #invtMM = ginv(t(M)%*%M)
  tXX = t(X)%*%X
  tUY = t(U)%*%Y
  tMX = t(M)%*%X

  #tUU = #rbind(cbind(tXX, t(tMX)),cbind(tMX, t(M)%*%M))
  #tUU.sqmat=sqrtmat.comp(tUU)
  tUU = ginv.largep(U,sqrtmat=TRUE,sqrtinvmat=TRUE)


  ## Interative Update

  #betaest =  matrix(0,1+2*V,length(lambda1)*length(lambda2)*length(tau)*length(alpha))
  lam1=rep(sort(lambda1,decreasing=TRUE), each=length(lambda2))
  lam2=rep(lambda2, length(lambda1))
  alpha = sort(alpha,decreasing=TRUE)

  myfunc<-function(j, k,gamma_new = rep(0,V+1),alpha_new = rep(0,V)){
    if(verbose==TRUE){print(paste("Lambda1=",lam1[j], "Lambda2=",lam2[j], "alpha=",alpha[k]))}

      iter=0
      err=1000
      allzero.count=0
      sigma2penalty=matrix(1,V,V);diag(sigma2penalty)<-0

      while( err>tol & iter<max.iter & allzero.count<4){

        alpha_old=alpha_new
        gamma_old = gamma_new
        beta_old = c(gamma_old,alpha_old)

        sigma1 = mean( (Y - U %*% gamma_old)^2)
        tmp = M - matrix(X,N,1) %*% matrix(alpha_old,1,V)
        Sigma2 = t(tmp)%*%tmp/N

        Omega=QUIC( Sigma2,rho=sigma2penalty*lam2[j],msg=0)#Inverse matrix of the covariance matrix of M

        Omega.sqrtmat=try(t(base::chol(Omega$X)),TRUE)
        if (is.matrix(Omega.sqrtmat)==FALSE){
          tmp.omega.1=base::chol(Omega$X,pivot=TRUE)
          Omega.sqrtmat=t(tmp.omega.1[,order(attr(tmp.omega.1, 'pivot'))])
        }
        #sqrtmat.comp(Omega$X)
        Omega.sqrtmat.inv=try(t(base::chol(Omega$W)),TRUE)#sqrtmat.comp(Omega$W)
        if (is.matrix(Omega.sqrtmat.inv)==FALSE){
          tmp.omega.2=base::chol(Omega$W,pivot=TRUE)
          Omega.sqrtmat.inv=t(tmp.omega.2[,order(attr(tmp.omega.2, 'pivot'))])
        }


        Asqmat = bdiag(1/sqrt(sigma1) * tUU$sqrtmat,  sqrt(as.numeric(tXX)) * Omega.sqrtmat)
        Asqmat.inv=bdiag(sqrt(sigma1) * tUU$sqrtinv,  1/sqrt(as.numeric(tXX)) * Omega.sqrtmat.inv)
        C = Asqmat.inv %*% rbind(tUY/sigma1, Omega$X%*%tMX)


        # fit = SGL(list(x=as.matrix(Asqmat), y=as.matrix(C)),
        #           index=grpgroup,
        #           lambdas=lam1[j],
        #           alpha=alphalist[j],standardize=FALSE)
        fit = oneDim_allowcov(list(x=as.matrix(Asqmat), y=as.matrix(C)),
                              index=c(1, rep(1:V +1,2)),
                              lambdas=lam1[j],
                              alpha=alpha[k],
                              group.penalty.factor = group.penalty.factor,
                              penalty.factor=penalty.factor)#,standardize=FALSE)

        fit$beta[abs(fit$beta)<threshold]<-0
        beta_new=fit$beta
        if (all(beta_new[-1]==0)){allzero.count=allzero.count+1}
        gamma_new = beta_new[1:(V+1)]#beta_new[c(1, (1:V)*2)]#
        alpha_new = beta_new[(1:V)+ V+1]#beta_new[c(1:V)*2+1]#beta_new[(1:V)+ V+1]#

 #       print(cbind(alpha_new,gamma_new[-1])[1:7,])
        err = sqrt(sum((beta_old[-1]-c(gamma_new[-1],alpha_new))^2))
        iter=iter+1
        if (verbose==TRUE){print(c(iter, err))}
      }
      ### compute BIC
      # zerolist=(c(gamma_new[1],alpha_new) ==0)
      # tmp = M - matrix(X,N,1) %*% matrix(alpha_new,1,V)
      # Sigma2 = t(tmp)%*%tmp/N
      # Omega=QUIC( Sigma2,rho=sigma2penalty*lam2[j],msg=0)
      # bic=N*log(sum(Y - cbind(X,M) %*% gamma_new)^2/N) + N*log(det(Omega$W)) +
      #   log(N)*(sum(1-zerolist))

      return(list(betahat=beta_new,Omegahat=Omega,sigma_y_sq=sigma1,#bic=bic,
                  alpha=alpha[k], lambda1=lam1[j],lambda2=lam2[j]))
  }

  zzz<-list()

  ## when the algorithm selects too many parameters, we stop there.
  for (k in 1:length(alpha)){
    zzz[[k]]<-list()
    j=0
    nonzeros=0
    while( j<length(lam1) & nonzeros< non.zeros.stop){
      j=j+1
      re<-c();
      gamma_init = rep(0,V+1)
      alpha_init = rep(0,V)
      zzz[[k]][[j]]<-NULL
      if (j>1){if(is.null(zzz[[k]][[j-1]])==FALSE){
        gamma_init = zzz[[k]][[j-1]]$betahat[1:(V+1)]
        alpha_init = zzz[[k]][[j-1]]$betahat[1:(V) + (V+1)]
      }}
      try(re<-myfunc(j=j,k=k,gamma_new=gamma_init,alpha_new=alpha_init))
      zzz[[k]][[j]]<-re
      nonzeros=sum(re$betahat!=0)
    }
  }

  betaest=do.call(cbind,lapply(zzz, function(x0){do.call(cbind,lapply(x0,function(xx){xx$betahat}))}))
  alphas=unlist(lapply(zzz, function(x0){unlist(lapply(x0,function(xx){xx$alpha}))}))
  sigma_y_sq=unlist(lapply(zzz, function(x0){unlist(lapply(x0,function(xx){xx$sigma_y_sq}))}))
  lam1s=unlist(lapply(zzz, function(x0){unlist(lapply(x0,function(xx){xx$lambda1}))}))
  lam2s=unlist(lapply(zzz, function(x0){unlist(lapply(x0,function(xx){xx$lambda2}))}))

  cest =betaest[1,]
  medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]
  nump=apply(betaest,2,function(x){sum(abs(x)>0)})

  if (Omega.out==FALSE){Omega=NULL;
  }else{Omegas=lapply(zzz, function(x0)lapply(x0,function(xx)xx$Omegahat$X));Omega=do.call(c, Omegas)}
  return(list(
    c = cest,
    hatb=betaest[(1:V)+1,]*Y.sd/M.sd,
    hata=betaest[(1:V)+V+1,]*M.sd/X.sd,
    medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]*Y.sd/X.sd,
    alpha=alphas,
    lambda1 = lam1s,
    lambda2= lam2s,
    nump=nump,
    Omega=Omega,
    sigma_y_sq=sigma_y_sq,
#    bic=bics,
    nmed=apply(as.matrix(medest), 2,function(x)sum(abs(x)>0))
  ))
}
