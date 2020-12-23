#' This function estimates the variance-covariance matrix of Y.
#' @param Y is a list 
#' @param t is the vector of distinct time points at which Y has been observed.
#' @param t.star is the list of distinct time points at which Y has been observed.
#' @param index is the list of indices corresponding to the positions of the time points at which the data is available
#' @details Functional Principal Component Analysis has been performed to get the variance covariance matrix for the vector obtained from Y.
#' @return sig is the variance-covariance matrix of Y.
#' @return sigi is the inverse of the variance-covariance matrix of Y.
#' @export
#' @author Mityl Biswas

est_sig <- function(Y, t, t.star, index)
{
  require(face)
  l <- length(unlist(t.star))
  m <- length(t)
  m.i <- unlist(lapply(Y, function(x) length(x)))
  n <- length(Y)
  seq <- seq(1:n)
  dat <- data.frame(argvals = unlist(t.star), y = unlist(Y))
  dat$subj <- unlist(sapply(seq, function(x) rep(x, m.i[x])))
  f <- face.sparse(dat, argvals.new = t)
  eval <- f$eigenvalues
  evec <- f$eigenfunctions
  f2 <- f$var.error.hat[1]
  if(length(eval) == 1)
  {  sig <- eval * evec %*% t(evec) + f2*diag(m) 
  }else
    sig <- evec%*%diag(eval)%*%t(evec) + f2*diag(m)
  #abc <- sig[unlist(index),unlist(index)]
  #prev line gave matrix with non 0 off diag. Can we fix this?
  sigi <- solve(sig)
  ABC <- matrix(0, l, l)
  ub <- 0
  for( i in 1:n)
  {
    lb <- ub + 1
    ub <- ub + length(unlist(index[i]))
    ABC[lb:ub, lb:ub] <- sigi[unlist(index[i]), unlist(index[i])]
  }
  #as.matrix(bdiag(list(matrix(seq(1:4),2), matrix(seq(1:9),3))))
  return(list(sig = sig, sigi = ABC))
}