#' Compute inverse, squareroot and inverse of the square root of the covariance
#'
#' @param mat One-dimensional predictor
#' @param thresh (default=10^{-20}) Threshold for eigen values
#' @return squareroot of matrix
#' @examples
#' set.seed(1234)
#' x=scale(matrix(100*20,100,20),center=TRUE, scale=FALSE)
#' ginv.largep(x)
#'
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords internal
#' @import MASS
#' @export

ginv.largep<-function(x.c,sqrtmat=TRUE, sqrtinvmat=TRUE, thresh=10^{-20}){
  #n=nrow(x)
  #x.c=scale(x, center=TRUE, scale=FALSE)
  xxt.inv= ginv( x.c %*% t(x.c))
  tmp = xxt.inv %*% x.c
  sqrt.mat=sqrt.invmat=NULL
  if (sqrtinvmat==TRUE){
    sqrt.mat=t(sqrtmat.comp(xxt.inv,thresh=thresh) %*% x.c) %*% x.c
  }
  if (sqrtinvmat==TRUE){
    sqrt.invmat=t(sqrtmat.comp(xxt.inv,thresh=thresh) %*% x.c) %*% xxt.inv %*% x.c
  }
  return(list(inv=t(tmp) %*% tmp, sqrtinv=sqrt.invmat, sqrtmat=sqrt.mat))
}
