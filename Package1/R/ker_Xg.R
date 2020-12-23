#' This function provides the gaussian kernel corresponding to X and E
#' @param X is a matrix. Each column of X corresponds to q many covariates of an individual
#' @param E is the matrix corresponding to the component of the kernel that does not involve X.
#' @details The kronecker product of the kernel for X and the kernel that does not depend on X.
#' @return K is the kernel matrix that would have been if data were observed at all the available time points.
#' @export
#' @author Mityl Biswas
ker_Xg <- function(X, E)
{
  require(kernlab)
  kern <- rbfdot(0.5)
  seq <- seq(1:ncol(X))
  # X.star <- outer(seq, seq, function(i, j) as.numeric(kern(X[,i], X[,j])))
  X1 <- sapply(seq, function(j) sapply(seq, function(i) kern(X[,i], X[,j])))
  K <- kronecker(X1, E)
  return(K)
}
