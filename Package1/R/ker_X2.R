#' This function provides the quadratic kernel corresponding to X and E
#' @param X is a matrix. Each column of X corresponds to q many covariates of an individual
#' @param E is the matrix corresponding to the component of the kernel that does not involve X.
#' @details The kronecker product of the kernel for X and the kernel that does not depend on X.
#' @return K is the kernel matrix that would have been if data were observed at all the available time points.
#' @export
#' @author Mityl Biswas
ker_X2 <- function(X, E)
{
  X.star <- (1 + t(X) %*% X)^2
  K <- kronecker(X.star, E)
  return(K)
}