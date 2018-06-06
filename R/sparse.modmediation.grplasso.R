#' Conduct sparse moderated mediation with group LASSO
#'
#' Fit a mediation model via penalized maximum likelihood and structural equation model.
#' The regularization path is computed for the lasso or elasticnet penalty at a grid of
#' values for the regularization parameter lambda. Currently, mediation analysis is developed based on gaussian assumption.
#'
#' Multiple Mediaton Model:
#' (1) M = Xa1 + XZ a2 + z a3 + e1
#' (2) Y = Xc1 + XZ c2  + Mb1 + MZ b2 + Z c3 + e2
#' And in the optimization, we do not regularize c', due to the assumption of partial mediation.
#' @param X One-dimensional predictor
#' @param M Multivariate mediator
#' @param Y Outcome
#' @param Z Scalar moderator
#' @param tol (default -10^(-10)) convergence criterion
#' @param max.iter (default=100) maximum iteration
#' @param lambda (default=log(1+(1:50)/125)) tuning parameter for L1 penalization
#' @param grpgroup (default=c(rep(1,3),rep( 1:V +1,5)))
#' @param penalty.factor (default=c(0,rep(sqrt(2),V))) give different weight of penalization for the 2V mediation paths.
#' @param X.cont (defult=TRUE) Related to standardization. If X is categorical data, standardiziagn won't apply.
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
#' fit=sparse.mediation.grplasso(X,M,Y,verbose=FALSE, lambda = log(1+(1:20)/50))
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords highdimensional mediation glmnet
#' @import parallel
#' @import MASS
#' @import gglasso
#' @export
sparse.modmediation.grplasso = function(X,M,Y,Z,tol=10^(-10),max.iter=100,
                                     lambda = log(1+(1:50)/125),X.cont=TRUE,
                                     grpgroup=c(rep(1,3), rep(1:V+1,5)),
                                     penalty.factor=c(0,rep(1,V)),
                                     threshold=0,
                                     verbose=FALSE){

  ## All variables has to be centered to zero and scaled to have variance 1.


  ## Center all values, and also make their scales to be 1. In this context, all coefficients will be dexribed in terms of correlation or partial correlations.
  N = nrow(M)
  V = ncol(M)

  Y = scale(Y,center=TRUE,scale=TRUE)
  if(X.cont==TRUE){X = matrix(scale(X,center=TRUE,scale=TRUE),N,1)}
  if(X.cont==FALSE){X = matrix(X,N,1)}

  M = scale(M, center=TRUE,scale=TRUE)
  Z = matrix(scale(Z, center=TRUE,scale=TRUE),N,1)

    #Y.mean=mean(Y)
  #X.mean=mean(X)
  #M.mean=apply(M,2,mean)

  #Y.sd=as.vector(sqrt(var(Y)))
  #X.sd=as.vector(sqrt(var(X)))
  #M.sd=sqrt(apply(M,2,var))

  ## Penalty Factor
  if (ncol(X)>1){stop("X has more than 1 colum. Stop.")}
  ## Initialization###
  ## OLS Estimation ###
  U = cbind(X, Z*X, Z, M, as.vector(Z)*M)
  tUU = t(U)%*%U
  tUU.sqmat=sqrtmat.comp(tUU)
  invtUU = ginv(tUU)
  #invtMM = ginv(t(M)%*%M)

  tXX = t(X)%*%X
  tZZ = t(Z)%*%Z
  tXZXZ = t(Z*X)%*%(Z*X)
  tUY = t(U)%*%Y
  tMX = t(M)%*%X
  tMXZ = t(M)%*% (Z*X)
  ## Interative Update

  betaest =  matrix(0,5*V+3,length(lambda))
  for (j in 1:length(lambda)){
    if (verbose==TRUE){print(paste("Lambda",lambda[j]))}
    c_new= rnorm(3)#invtUU %*% tUY
    b_new=rnorm(2*V)
    a1_new=a2_new=a3_new = rnorm(V)#t(ginv(t(X)%*%X)%*%t(X)%*%M)

    iter=0
    err=1000
    while( err>tol & iter<max.iter){

      c_old=c_new
      b_old=b_new
      a1_old=a1_new
      a2_old=a2_new
      a3_old=a3_new

      beta_old = c(c_old,b_old,a1_old,a2_old,a3_old )

      sigma1 = mean((Y - cbind(X, Z*X, Z) %*% c_old - cbind(M, as.vector(Z)*M)%*%b_old )^2)
      tmp = M - cbind(X, Z*X, Z)%*%t(cbind(a1_old,a1_old,a1_old))
      Sigma2 = t(tmp)%*%tmp/N
      Sigma2.sqmat=sqrtmat.comp(Sigma2)
      Sigma2.sqrt.inv=ginv(Sigma2.sqmat)
      Sigma2.inv=ginv(Sigma2)


      #A = bdiag(1/sigma1 * tUU,bdiag(tXX, tXZXZ, tZZ) %x% Sigma2.inv)
    #  A[1:(3+2*V),1:(3+2*V)]=1/sigma1 * tUU
    #  A[(1:V)+(3+2*V),1:V+(3+2*V)]=  as.numeric(tXX) * Sigma2.inv
    #  A[(1:V)+(3+3*V),1:V+(3+3*V)]=  as.numeric(tXZXZ) * Sigma2.inv
    #  A[(1:V)+(3+4*V),1:V+(3+4*V)]=  as.numeric(tZZ) * Sigma2.inv
      sqmatA=bdiag(1/sqrt(sigma1) * tUU.sqmat,bdiag(sqrt(tXX), sqrt(tXZXZ), sqrt(tZZ)) %x% Sigma2.sqrt.inv)
      #print(sqmatA)
      #sqmatA = A;sqmatA[1:(1+V),1:(1+V)]=1/sqrt(sigma1) * tUU.sqmat
      #sqmatA[(1+V)+ 1:V,(1+V)+ 1:V]=  sqrt(as.numeric(tXX)) * Sigma2.sqrt.inv
      #C = ginv(sqmatA) %*% rbind(tUY/sigma1, Sigma2.inv%*%tMX, )
      C = ginv(as.matrix(sqmatA)) %*% c(tUY/sigma1, as.vector(Sigma2.inv%*%t(M)%*%cbind(X, Z*X, Z)))
      if(is.null(penalty.factor)==TRUE){
        #fit = glmnet(sqmatA, C,lambda=lambda[j],alpha=alpha)
        fit=gglasso(x=scale(sqmatA)[,order(grpgroup)], y=scale(C),lambda=lambda[j],group=grpgroup[order(grpgroup)])
      }else{
        #fit = glmnet(sqmatA, C,lambda=lambda[j],penalty.factor=penalty.factor,alpha=alpha)
        fit=gglasso(x=scale(sqmatA)[,order(grpgroup)], y=scale(C),lambda=lambda[j],group=grpgroup[order(grpgroup)],pf = penalty.factor)
      }

      beta_new = as.vector(coef(fit))[-1]
#print(beta_new)
      ## use thresholds as well: since all variables are standardized, coefficients less than 0.001 does not have any meaning.
      if (threshold>0){
        if (sum(abs(beta_new)<threshold)>0){
          beta_new[abs(beta_new)<threshold]<-0
        }
      }
#print(beta_new)
      c_new=beta_new[1:3]
      b_new=beta_new[ c((1:V)*5-1 ,(1:V)*5)]
      a1_new=beta_new[ (1:V)*5+1 ]
      a2_new=beta_new[ (1:V)*5+2 ]
      a3_new=beta_new[ (1:V)*5+3 ]

      beta_new<-beta_new[c(1:3,(1:V)*5-1 ,(1:V)*5,(1:V)*5+1,(1:V)*5+2,(1:V)*5+3 )]
      #print(beta_new)
      err = sum((beta_old[-c(1:3)]-beta_new[-c(1:3)])^2)
      iter=iter+1
      if (verbose==TRUE){print(c(iter, err))}
      #print(beta_new[1:15])
    }
    betaest[,j]=beta_new
  }
  cest =betaest[1:3,]
  medest = betaest[(1:V)+3,]*betaest[(1:V)+2*V+3,]
  modest=  betaest[(1:V)+3+V,]*betaest[(1:V)+3*V+3,]
  nump=apply(betaest,2,function(x){sum(abs(x)>0)})


  return(list(
    c = cest,
    hatb1=betaest[(1:V)+3,],
    hatb2=betaest[(1:V)+V+3,],
    hata1=betaest[(1:V)+2*V+3,],
    hata2=betaest[(1:V)+3*V+3,],
    hata3=betaest[(1:V)+4*V+3,],
    medest = medest,modest=modest,
    lambda = lambda,
    nump=nump
  ))
}

