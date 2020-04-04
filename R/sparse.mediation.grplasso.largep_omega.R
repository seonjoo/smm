#' Conduct sparse mediation for large p ( p > n) with grouplasso penalty and fast computation of inverse matrix
#'
#' Fit a mediation model via penalized maximum likelihood and structural equation model.
#' The regularization path is computed for the grouplasso penalty at a grid of
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
#' V=300
#' set.seed(1234)
#' a = c(rep(1,3),rep(0,V-3))*5;b<-a
#' X = rnorm(N)
#' M =  X %*% t(a)+ matrix(rnorm(N*V),N,V)
#' Y =  X + M %*% b + rnorm(N)
#' fit=sparse.mediation.grplasso.largep_omega(X,M,Y,verbose=FALSE)
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation L1penalization
#' @import parallel
#' @import MASS
#' @import glmnet
#' @import QUIC
#' @import Matrix
#' @importFrom stats var predict
#' @export
sparse.mediation.grplasso.largep_omega = function(X,M,Y,
                                                  tol=10^(-10),
                                                  max.iter=100,
                                                  lambda1 =exp(-8:0),
                                                  lambda2 = exp(seq(0,0.8*log(ncol(M)),length=4)),
                                                  grpgroup=c(1, rep(1:(ncol(M))+1,2)),
                                                  penalty.factor=c(0,rep(1,ncol(M))),
                                                  verbose=FALSE,
                                                  Omega.out=FALSE){
  ## Center all values, and also make their scales to be 1. In this context, all coefficients will be dexribed in terms of correlation or partial correlations.
  N = nrow(M)
  V = ncol(M)
  #Y.mean=mean(Y)
  #X.mean=mean(X)
  #M.mean=apply(M,2,mean)
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

  #invtMM = ginv(t(M)%*%M)
  tXX = t(X)%*%X
  tUY = t(U)%*%Y
  tMX = t(M)%*%X

  #tUU = #rbind(cbind(tXX, t(tMX)),cbind(tMX, t(M)%*%M))
  #tUU.sqmat=sqrtmat.comp(tUU)
  tUU = ginv.largep(U,sqrtmat=TRUE,sqrtinvmat=TRUE)


  ## Interative Update
  lam1=rep(lambda1, each=length(lambda2))
  lam2=rep(lambda2, length(lambda1))
#  for (j in 1:length(lam1)){
  myfunc<-function(j){
    if(verbose==TRUE){print(paste("Lambda1=",lam1[j], "Lambda2=",lam2[j]))}
      gamma_new = rep(0,V+1)#tUU$inv %*% tUY
      alpha_new = rep(0,V)#t(ginv(tXX)%*%t(X)%*%M)

      iter=0
      err=1000
      while( err>tol & iter<max.iter){

        alpha_old=alpha_new
        gamma_old = gamma_new
        beta_old = c(gamma_old,alpha_old)

        sigma1 = mean( (Y - U %*% gamma_old)^2)
        tmp = M - matrix(X,N,1) %*% matrix(alpha_old,1,V)
        Sigma2 = t(tmp)%*%tmp/N

        Omega=QUIC( Sigma2,rho=lam2[j],msg=0)#Inverse matrix of the covariance matrix of M

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


        sqmatA = bdiag(1/sqrt(sigma1) * tUU$sqrtmat,  sqrt(as.numeric(tXX)) * Omega.sqrtmat)
        sqmatA.inv=bdiag(sqrt(sigma1) * tUU$sqrtinv,  1/sqrt(as.numeric(tXX)) * Omega.sqrtmat.inv)
        C = sqmatA.inv %*% rbind(tUY/sigma1, Omega$X%*%tMX)

        if(is.null(penalty.factor)==TRUE){
          #fit = glmnet(sqmatA, C,lambda=lambda[j],alpha=alpha)
          fit=gglasso(x=scale(sqmatA)[,order(grpgroup)],
                      y=scale(C),
                      lambda=lam1[j],
                      group=grpgroup[order(grpgroup)])
        }else{
          #fit = glmnet(sqmatA, C,lambda=lambda[j],penalty.factor=penalty.factor,alpha=alpha)
          fit=gglasso(x=scale(sqmatA)[,order(grpgroup)],
                      y=scale(C),
                      lambda=lam1[j],
                      group=grpgroup[order(grpgroup)],
                      pf = penalty.factor)
        }
        beta_new=coef(fit)[-1]
        gamma_new = beta_new[c(1, (1:V)*2)]#beta_new[1:(V+1)]
        alpha_new = beta_new[c(1:V)*2+1]#beta_new[(1:V)+ V+1]
        err = sum((beta_old[-1]-c(gamma_new[-1],alpha_new))^2)
        iter=iter+1
        if (verbose==TRUE){print(c(iter, err))}
      }
      return(list(betahat=beta_new[c(1, (1:V)*2, c(1:V)*2+1)],Omegahat=Omega))
  }
  zzz=lapply(1:length(lam1), function(xxx){
    re<-c();try(re<-myfunc(xxx));return(re)})

  betaest=do.call(cbind,lapply(zzz, function(x)x$betahat))
  #     betaest[,j]=beta_new
 #   }

  cest =betaest[1,]
  medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]
  nump=apply(as.matrix(betaest),2,function(x){sum(abs(x)>0)})

  if(Omega.out==FALSE){Omegas=NULL}
  if (Omega.out==TRUE){Omegas=lapply(zzz, function(x)x$Omegahat)}
  return(list(
    c = cest,
    hatb=betaest[(1:V)+1,]*Y.sd/M.sd,
    hata=betaest[(1:V)+V+1,]*M.sd/X.sd,
    medest = betaest[(1:V)+1,]*betaest[(1:V)+V+1,]*Y.sd/X.sd,
#    alpha=alphalist,
#    tau=taulist,
    lambda1 = lam1,
    lambda2 = lam2,
    nump=nump,
    Omega=Omegas
  ))
}
