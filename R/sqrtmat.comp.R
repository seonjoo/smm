#' Compute squareroot of large covariance matrix
#'
#' @param mat Covariance matrix
#' @param thresh (default: 10^(-20)) for matrix decomposition
#' @param K Dimension of the matrix. Default is set as ncol(mat)-1
#' @return squareroot of matrix
#' @examples
#' set.seed(1234)
#' x=matrix(100*20,100,20)
#' covmat=cov(x)
#' sqrtmat.comp(covmat)
#'
#' @author Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
#' @references TBA
#' @keywords internal
#' @import MASS
#' @import rsvd
#' @export

sqrtmat.comp<-function(mat,thresh=10^(-20),K=NULL){
  if(is.null(K)){K=ncol(mat)}
  if (ncol(mat)>200){
    eigenmat=rsvd(mat, k=K)
  }else{eigenmat=svd(mat, nv=K, nu=K)}
  #  print(paste('Dimension of mat:', dim(mat)))
  ncomp=sum(eigenmat$d>thresh)
  #print(ncomp)
  # print(eigenmat$d)
  if (ncomp<2){
    sqmat=as.matrix(eigenmat$v[,1]) %*% sqrt(eigenmat$d[1]) %*% t(as.matrix(eigenmat$v[,1]))
  }else{sqmat=eigenmat$v[,1:ncomp] %*% diag(sqrt(eigenmat$d[1:ncomp])) %*% t(eigenmat$v[,1:ncomp])
  }


  return(sqmat)
}
