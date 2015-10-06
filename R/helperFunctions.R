#' Check if a matrix is computationally positive-definite
#' @param M matrix of interest
#' @return TRUE if matrix is computationally positive-definite
isComputationallyPD <- function(M){
 eig <- eigen(M, only.values=TRUE, symmetric=TRUE)
 if(min(eig$values)>0){
    # All eigenvalues >0 
    if(rcond(M)>.Machine$double.eps^2){
      # Has a condition number greater than machine precision
      PD <- TRUE
    } else {
      # Fails if condition number too small
      PD <- FALSE
    }
 } else {
   # Fails on basis of eigenvalues
   PD <- FALSE
 }
 # Return PD status
 return(PD)
}