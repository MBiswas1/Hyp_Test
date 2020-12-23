#' This function provides the test statistic
#' @param Y.sigi is a vector comprising the centered observations stacked by individuals.
#' @param K is the kernel matrix.
#' @details The test statistic has not been centered.
#' @return T is the desired test statistic.
#' @export
#' @author Mityl Biswas
test_stat <- function(Y.sigi, K)
{
  T <- t(Y.sigi) %*% K %*% Y.sigi
  return(T)
}