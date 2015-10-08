#' Test for a perturbation at a specific location
#' 
#' @param Y matrix or data frame (rows=subjects, columns=variables)
#' @param Omega precision matrix
#' @param Sigma covariance matrix
#' @param perturb.vec boolean vector indicating which variables to test for perturbations
#' @param return.value value to return (\code{'stat'} or \code{'pval'})
#' @return test statistic for a perturbation at perturb.vec=TRUE locations only
#' 
#' @references
#' Griffin, P.J, Johnson, W.E, and Kolaczyk, E.K. Detection of multiple 
#' perturbations in multi-omics biological networks (under review)
#' 
#' @noRd
perturbTest <- function(Y, Omega, Sigma, perturb.vec, return.value='stat'){
  perturb.vec <- as.logical(perturb.vec)
  vec. <- which(!perturb.vec)
  j <- which(perturb.vec)
  
  n <- nrow(Y)
  y <- colMeans(Y)
  z <- Omega %*% y
  R <- chol(Sigma)
  l.null <- n*crossprod(R %*% z)
  
  if(sum(vec.)>0){
    R11 <- chol(Sigma[j, j, drop=FALSE])
    half <- solve(t(R11), Sigma[j, vec., drop=FALSE])
    inner <- Sigma[vec., vec., drop=FALSE] - crossprod(half)
    R.alt <- chol(inner)
    l.alt  <- n*crossprod(R.alt %*% z[vec.])
  } else {
    l.alt <- 0
  }
  
  # Return test statistic or p-value as requested
  if( return.value == 'stat'){
    retval <- l.null - l.alt
  } else if( return.value == 'pval'){
    retval <- pchisq(l.null - l.alt, df=sum(perturb.vec != 0), lower.tail=FALSE)
  } else {
    print('Error: return.value must be either stat or pval') 
  }
  return(retval)
}


#' Test for perturbations at a series of locations
#' 
#' Given a matrix of possible perturbation locations, perform the single-target
#' (non-sequential) testing procedure described by Griffin et al. Return either
#' test statistics or p-values.
#' 
#' @param Y matrix or data frame (rows=subjects, columns=variables)
#' @param Omega precision matrix
#' @param Sigma covariance matrix
#' @param perturb.mat matrix indicating sets of which variables to test for perturbations
#' @param return.value value to return (\code{'stat'} or \code{'pval'})
#' @return test statistics for perturbations at each of perturb.mat columns
#' 
#' @references
#' Griffin, P.J., Johnson, W.E., and Kolaczyk, E.K. Detection of multiple 
#' perturbations in multi-omics biological networks (under review)
#' 
#' @importFrom mvtnorm rmvnorm
#' 
#' @examples
#' library(mvtnorm)
#' Omega <- matrix(0, nrow=6, ncol=6)
#' Omega[1:4,1:4] <- -1
#' diag(Omega) <- 6
#' Sigma <- solve(Omega)
#' mu <- c(1,1,0,0,0,0)
#' Y <- rmvnorm(20, mean=Sigma %*% mu, sigma=Sigma)
#' perturb.mat <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
#'                         FALSE, FALSE, TRUE, TRUE, FALSE, FALSE), ncol=2)
#' perturbTests(Y, Omega, Sigma, perturb.mat, return.value='stat')
#' 
#' @export
perturbTests <- function(Y, Omega, Sigma, perturb.mat, return.value='stat'){  
  retval <- apply(perturb.mat, 2, perturbTest, Y=Y, Omega=Omega, Sigma=Sigma, return.value=return.value)
  return(retval)
}


#' Sequentially test for perturbations at a series of locations
#' 
#' Given a matrix of possible perturbation locations, perform the sequential
#' testing procedure described by Griffin et al. Return either test statistics
#' or p-values.
#' 
#' @param Y matrix or data frame (rows=subjects, columns=variables)
#' @param Omega precision matrix
#' @param Sigma covariance matrix
#' @param perturb.mat matrix indicating sets of which variables to test for perturbations
#' @param return.value return statistics (\code{'stat'}) or p-values (\code{'pval'})
#' @return test statistics for perturbations at each of perturb.mat columns
#' 
#' @references
#' Griffin, P.J., Johnson, W.E., and Kolaczyk, E.K. Detection of multiple 
#' perturbations in multi-omics biological networks (under review)
#' 
#' @examples
#' library(mvtnorm)
#' Omega <- matrix(0, nrow=6, ncol=6)
#' Omega[1:4,1:4] <- -1
#' diag(Omega) <- 6
#' Sigma <- solve(Omega)
#' mu <- c(1,1,0,0,0,0)
#' Y <- rmvnorm(20, mean=Sigma %*% mu, sigma=Sigma)
#' perturb.mat <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
#'                         FALSE, FALSE, TRUE, TRUE, FALSE, FALSE), ncol=2)
#' perturbTests(Y, Omega, Sigma, perturb.mat, return.value='stat') # initial ranking            
#' seqPerturbTests(Y, Omega, Sigma, perturb.mat, return.value='stat') # sequential adjusted
#' 
#' @export
seqPerturbTests <- function(Y, Omega, Sigma, perturb.mat, return.value='stat'){
  # First test is conditional on nothing
  lrt.init <- perturbTests(Y, Omega, Sigma, perturb.mat, 'stat')
  col.null <- which.max(lrt.init)
  col.stat <- max(lrt.init)
  test.mat <- perturb.mat
  
  # For remaining spaces, conditional tests
  for(i in 1:(ncol(perturb.mat)-1) ) {
    found.locs <- rowSums(perturb.mat[,col.null,drop=FALSE])>0
    test.mat[found.locs,] <- TRUE
    
    null.test <- c(perturbTest(Y, Omega, Sigma, found.locs, 'stat'))
    alt.test <- perturbTests(Y, Omega, Sigma, test.mat, 'stat')
    lrt.cond <- (alt.test - null.test)
    
    col.null <- c(col.null, which.max(lrt.cond))
    col.stat <- c(col.stat, max(lrt.cond))
  }
  stats <- rep(NA, length(col.null))
  stats[col.null] <- col.stat
  
  # Convert to p-values if requested
  if(return.value=='pval'){
    retval <- pchisq(stats, df=colSums(perturb.mat), lower.tail=FALSE)
  } else if (return.value=='stat') {
    retval <- stats 
  } else {
    print('Error: return.value must be either stat or pval') 
  }
  
  return(retval)
}

#' Get node-wise test statistics 
#' 
#' Given a vector \code{id} indicating node membership, test for perturbations 
#' at each node as described by Griffin et al.  This function performs either 
#' the single-target test (\code{sequential=FALSE}, default) or the 
#' multi-target adjusted procedure  (\code{sequential=TRUE})
#' 
#' @param Y matrix or data frame (rows=subjects, columns=variables)
#' @param Omega precision matrix
#' @param Sigma covariance matrix
#' @param id vector mapping variables to nodes
#' @param sequential boolean, whether or not to perform sequential adjustments
#' @param return.value return statistics (\code{'stat'}) or p-values (\code{'pval'})
#' @return test statistics or p-values for perturbations at each individual node
#' 
#' @references
#' Griffin, P.J., Johnson, W.E., and Kolaczyk, E.K. Detection of multiple 
#' perturbations in multi-omics biological networks (under review)
#' 
#' @examples
#' library(mvtnorm)
#' Omega <- matrix(0, nrow=6, ncol=6)
#' Omega[1:4,1:4] <- -1
#' diag(Omega) <- 6
#' Sigma <- solve(Omega)
#' mu <- c(1,1,0,0,0,0)
#' Y <- rmvnorm(20, mean=Sigma %*% mu, sigma=Sigma)
#' id <- c(1,1,2,2,3,3)
#' perturbNodeTests(Y, Omega, Sigma, id, sequential=FALSE) # unadjusted           
#' perturbNodeTests(Y, Omega, Sigma, id, sequential=TRUE) # sequential
#' 
#' @export
perturbNodeTests <- function(Y, Omega, Sigma, id, sequential=FALSE, 
                              return.value='stat'){
  p <- length(unique(id))
  node.perturb <- 1:p
  perturb.mat <- sapply(id, function(i)(id==i))
  dup.cols <- duplicated(t(perturb.mat))
  perturb.mat <- perturb.mat[,!dup.cols] # one column per node to test.
  
  if(!sequential){
    retval <- perturbTests(Y=Y, Omega=Omega, Sigma=Sigma, perturb.mat=perturb.mat, 
                     return.value=return.value)
  } else {
    retval <- seqPerturbTests(Y=Y, Omega=Omega, Sigma=Sigma, perturb.mat=perturb.mat,
                         return.value=return.value)
  }
  return(retval)
}
