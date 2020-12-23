#' This function provides the kernel for the observed points
#' @param index.star is the vector of indices for the time points at which the data was observed. 
#' @param K is the matrix corresponding to the kernel that would have been if the data were dense.
#' @details We obtain the kernel matrix at the required indices.
#' @return K.star is the kernel matrix
#' @export
#' @author Mityl Biswas
kern_miss <- function(index.star, K)
{
  K.star <- K[index.star,index.star]
  return(K.star)
}