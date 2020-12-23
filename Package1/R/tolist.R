#' This function converts a vector to a list.
#' @param Y is the vector that needs to be converted to a list.
#' @param m.i is a vector comprising the lengths of each element of the resulting list.
#' @details The inverse operation of the unlist function.
#' @return Y.l is Y as a list.
#' @export
#' @author Mityl Biswas
tolist <- function(Y, m.i)
{
  n <- length(m.i)
  c <- 1
  Y.l <- NULL
  for(i in 1:n)
  {
    Y.l[[i]] <- Y[c : (c + m.i[i] - 1)]
    c <- c + m.i[i]
  }
  return(Y.l)
}