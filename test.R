fit<-sparse.mediation.grplasso(x,m,y,tol=0.0001,max.iter=50,
                               lambda1 =exp(seq(-2,-5,length=10)),
                               lambda2=1,alpha=c(0.2,0.5),verbose=TRUE)

library(smm)
rm(list=ls())
N=200
V=100
set.seed(1234)
a = c(rep(1,3),rep(0,V-3));
b1=b2=rep(0,V)
b1[1]<-1
b2[1]<-1
x = rnorm(N)
m =  x %*% t(a)+ matrix(rnorm(N*V),N,V)
y =  x +m %*% b1 + (x*m) %*% b2 +  rnorm(N)
X=x;M=m;Y=y
fit=sparse.txtmedint.sgrlasso.largep_omega(x,m,y,alpha=c(0.55,0.75,.95),
                                           lambda1=exp(seq(-0.5, -3, length = 10)),
                                           non.zeros.stop=50,
                                           group.penalty.factor=c(0,rep(1,V)),
                                           penalty.factor=c(0,rep(1, 3*V)),
                                           verbose=TRUE,tol=0.0001,
                                           max.iter=50)
fit$c
fit

fit=sparse.mediation.sgrlasso.largep_omega(x,m,y,alpha=c(0.55,0.75,.95),
                                           lambda1=exp(seq(0.5, -3, length = 10)),
                                           non.zeros.stop=100,
                                           group.penalty.factor=c(0,rep(1,V)),
                                           penalty.factor=c(0,rep(1, 2*V)),
                                           verbose=TRUE,tol=0.0001,
                                           max.iter=50)
fit$c
fit
###################################
N=200
V=200
set.seed(1234)
a = c(rep(0.5,3),rep(0,V-3));
b1=b2=rep(0,V)
b1[1]<-1
x = rnorm(N)
m =  x %*% t(a)+ matrix(rnorm(N*V),N,V)
y =  x +m %*% b1  +  rnorm(N)

fit.cv = cv.sparse.mediation.sgrplasso(x,m,y,
                                       K=4,
                                       multicore=4,
                                       lambda1=exp(seq(0.5, -3, length = 10)),
                                       lambda2=exp(-1),
                                       alpha=c(0.55,0.75,.95),
                                       non.zeros.stop=100,
                                       group.penalty.factor=c(0,rep(1,V)),
                                       penalty.factor=c(0,rep(1, 2*V)),
                                       tol=0.0001,
                                       max.iter=50)

fit=sparse.mediation.grplasso(x,m,y,
                              lambda1=fit.cv$cv.lambda1,
                              lambda2=fit.cv$cv.lambda2,
                              alpha=fit.cv$cv.alpha,
                              group.penalty.factor=c(0,rep(1,V)),
                              penalty.factor=c(0,rep(1, 2*V)),
                              verbose=TRUE,tol=0.0001,
                              max.iter=50)
fit
