qrsh -l mem=10G,time=1:: -pe smp 8
/nfs/apps/R/3.6.0/bin/R

#.libPaths('/ifs/scratch/msph/LeeLab/softwares/R/hpc')
#library(QUIC)
#library(Matrix)
rm(list=ls())
library(smm)
library(MASS)
library(ppcor)
################################################################################

datget<- function(seednum=213456,N=200,V=800,
                 lambda1=exp(seq(0,-6, length=15)),
                 lambda2=exp(c(-1,1)),
                 alpha=c(0.5,0.75,0.95)){
  cmat=matrix(0,V,V);
  a=c(0.5,0.4,0.3)
  totef=0.5
  b=a*(totef+sqrt((1-totef^2)*(1-a^2)))
  cmat[1,2]<-totef
  cmat[1,1:3+2]<-a
  cmat[2,1:3+2]<-b
  cmat=cmat+t(cmat)
  diag(cmat)<-1

  set.seed(seednum)
  print(seednum)
  v=ncol(cmat)-2
  tmpmat = mvrnorm(n=N, mu=rep(0,v+2),Sigma=cmat)
  x=as.matrix(tmpmat[,1])
  m=as.matrix(tmpmat[,1:v + 2])
  y= as.matrix(tmpmat[,2])
### minimum 200 and V/4
  return(c(pcor.test(m[,1],y,x)$estimate,pcor.test(m[,2],y,x)$estimate,pcor.test(m[,3],y,x)$estimate))
}

aaa=do.call(cbind,lapply(1:100,function(x)datget(seednum=x)))
round(apply(aaa,1,mean),2)
# library(mediation)
# tmpdat=data.frame(y=y,m=m[,1:3], x=x)
# fit1=lm(m.3 ~ x,tmpdat)
# fit2=lm(y~ x + m.1 + m.2+ m.3 , tmpdat)
# med.fit=mediate(fit1, fit2, mediator='m.3',treat='x')
# summary(med.fit)

system.time(
fit.cv<-cv.sparse.mediation.grplasso(x,m,y,tol=tol,K=K,max.iter=max.iter,
                                         lambda1 = lambda1,lambda2=lambda2,
                                         alpha=alpha,multicore=multicore,
                                         non.zeros.stop=max(150,V/2)))

fit.cv$cv.lambda1
fit.cv$cv.lambda2
fit.cv$cv.alpha

system.time(
    fit<-sparse.mediation.grplasso(x,m,y,tol=tol,max.iter=max.iter,
                                                    lambda1 =fit.cv$cv.lambda1,
                                                    lambda2=fit.cv$cv.lambda2,
                                                    alpha=fit.cv$cv.alpha,
                                                    non.zeros.stop=max(150,V/2),verbose=TRUE))


fit
fit$med[1:3]
