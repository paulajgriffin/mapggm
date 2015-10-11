#' Get Frobenius norms of submatrices with optional weights
#' 
#' Given a square matrix \code{M} and a vector \code{id} that distinguishes sections of \code{M}, 
#' return a matrix with the Frobenius norm of the specified submatrices.  
#' Optionally, provide a square weight matrix \code{W} with number of rows and columns
#' equal to the unique entries in id that multiplies these norms.
#' 
#' @param M square matrix of interest
#' @param id vector grouping elements of M 
#' @param W optional weights (square matrix, dimensions equal to unique id
#' @return matrix of Frobenius norms of submatrices
#' 
#' @export
blockNorms <- function(M, id, W=NULL){
  # Weights 
  nodes <- unique(id)
  n.nodes <- length(unique(id))
  if (is.null(W)) {
    W <- matrix(1, n.nodes, n.nodes) 
  } else if( any(dim(W) != rep(length(nodes),2)) ) {
   stop('Error: incompatible W and id') 
  }
  rownames(W) <- colnames(W) <- nodes
  
  fnorms <- matrix(0, ncol=ncol(W), nrow=nrow(W), dimnames=list(nodes, nodes) )
  for (a in nodes){
   for(b in nodes){
     fnorms[which(nodes==a),which(nodes==b)] <- W[which(nodes==a), 
                                                  which(nodes==b), drop=FALSE] *
       norm(M[which(id==a), which(id==b), drop=FALSE], type='F')
   }
  }
  # Return the matrix of norms
  return(fnorms)
}

#' Get BIC of a given network configuration
#' 
#' In a zero-mean Gaussian graphical model, calculate the Bayesian information
#' criterion (BIC) of a model with the proposed precision matrix Omega.
#' 
#' @param S sample covariance matrix
#' @param n number of samples
#' @param Omega estimated precision
#' @param id vector grouping elements of S, Omega
#' @return Bayesian information criterion for model implied by Omega
#' 
#' @export
getBIC <- function(S, n, Omega, id){  
  nodes <- unique(id)
  n.nodes <- length(unique(id))
  
  blocks <- matrix(0, n.nodes, n.nodes)
  for (a in nodes){
   for (b in nodes[1:length(nodes)<=(which(nodes==a))]){
     blocks[which(nodes==a),which(nodes==b)] <- 
       ifelse(max(abs(Omega[which(id==a), which(id==b)]))!=0, 
              sum(id==a)*sum(id==b), 0)
   }
  }
  blockDims <- sum(blocks)
  
  BIC <- n*(sum(diag(S %*% Omega)) - log(det(Omega))) +  blockDims * log(n)
  return(BIC)
}


#' Get extended BIC of a given network configuration
#' 
#' In a zero-mean Gaussian graphical model, calculate the extended Bayesian 
#' information criterion (EBIC; Chen and Chen, 2008) of a model with the 
#' proposed precision matrix Omega.  The parameter gamma controls the weight
#' given towards the component of EBIC related to model size.
#' 
#' @param S sample covariance matrix
#' @param n number of samples
#' @param Omega estimated precision
#' @param id vector grouping elements of S, Omega
#' @param gamma EBIC parameter (default 0.5)
#' @return extended Bayesian information criterion for model implied by Omega
#' 
#' @references
#' Chen, J. and Chen, Z. (2008). Extended Bayesian information criteria for 
#' model selection with large model spaces. Biometrika 95, 759-771.
#' 
#' @export
getEBIC <- function(S, n, Omega, id, gamma=0.5){  
  nodes <- unique(id)
  n.nodes <- length(unique(id))
  
  blocks <- matrix(0, n.nodes, n.nodes)
  for (a in nodes){
   for (b in nodes[1:length(nodes)<=(which(nodes==a))]){
     blocks[which(nodes==a),which(nodes==b)] <- 
       ifelse(max(abs(Omega[which(id==a), which(id==b)]))!=0, 
              sum(id==a)*sum(id==b), 0)
   }
  }
  blockDims <- sum(blocks)
  lid <- length(id)
  p <- lid*(lid-1)/2
  k <- 0:p
  x <- -gamma*lgamma(k+1)-gamma*lgamma(p-k+1)
  mx <- max(x)
  ldenom <- log(sum(sort(exp(x - mx)))) + gamma*lgamma(p+1) + mx
  lnum <- gamma*lchoose(p,k)
  lprobs <- lnum-ldenom
  
  EBIC <- n*(sum(diag(S %*% Omega)) - log(det(Omega)))  +  blockDims * log(n) + 
    2*gamma*lprobs[blockDims]
  
  return(EBIC)
}

#' Get AIC of a given network configuration
#' 
#' In a zero-mean Gaussian graphical model, calculate the Akaike information
#' criterion (AIC) of a model with the proposed precision matrix Omega.
#' 
#' @param S sample covariance matrix
#' @param n number of samples
#' @param Omega estimated precision
#' @param id vector grouping elements of S, Omega
#' @return Akaike information criterion for model implied by Omega
#' 
#' @export
getAIC <- function(S, n, Omega, id){  
  nodes <- unique(id)
  n.nodes <- length(unique(id))
  
  blocks <- matrix(0, n.nodes, n.nodes)
  for (a in nodes){
   for (b in nodes[1:length(nodes)<=(which(nodes==a))]){
     blocks[which(nodes==a),which(nodes==b)] <- 
       ifelse(max(abs(Omega[which(id==a), which(id==b)]))!=0, 
              sum(id==a)*sum(id==b), 0)
   }
  }
  blockDims <- sum(blocks)
  AIC <- n*(sum(diag(S %*% Omega)) -  log(det(Omega)) ) +  blockDims * 2
  return(AIC)
}

#' Get Hamming distance of a given network configuration (truth required)
#' 
#' Calculate the Hamming distance of a given network construction (as provided
#' by Omega), compared to a true network.  
#' 
#' @param Omega estimated precision
#' @param Theta true network graph of joint nodes
#' @param id vector grouping elements of Omega
#' @return Hamming distance between Theta and graph implied by Omega
#' 
#' @export
getHammingDist <- function(Omega, Theta, id){
  # Overall node graph
  p <- length(unique(id))
  Theta.combo <- Theta[!duplicated(id), !duplicated(id)]
  ans.combo <- Matrix(ifelse(blockNorms(M=Omega, W=matrix(1,p,p ), 
                                        id=id)!=0,1,0), sparse=TRUE)
  diag(ans.combo) <- 0
  
  hamming <- sum(abs(ans.combo-Theta.combo))/2
  return(hamming)
}


#' Get value of primal function (called by maLasso1)
#' 
#' 
#' @param Omega current estimated precision
#' @param S sample covariance matrix
#' @param W weight matrix
#' @param id vector grouping elements of S, Omega
#' @param lambda penalty parameter
#' @return value of primal function for optimization problem
#' @noRd
getPrimal <- function(Omega, S, W, id, lambda){
  primal <- sum(diag(S %*% Omega)) - log(det(Omega)) + 
    lambda*sum(blockNorms(M=Omega, W=W, id=id))
  
  return(primal)
}

#' Get value of dual function (called by maLasso1)
#' 
#' @param Sigma current estimated covariance
#' @return value of dual function for optimization problem
#' @noRd
getDual <- function(Sigma) {
 dual <- ncol(Sigma) + log(det(Sigma))
 return(dual)
}

#' Update covariance/precision matrix estimates
#' 
#' @param t.step step size currently in use
#' @param id vector of node identifiers
#' @param S sample covariance
#' @param Omega.tmp current estimate of precision matrix
#' @param Sigma.tmp current estimate of covariance matrix
#' @param lambda penalty parameter
#' @param W weight matrix
#' @return list of updated covariance and precision estimates
#' @noRd
updateNodes <- function(t.step, id, S, Omega.tmp, Sigma.tmp, lambda, W ){
  
  # Blank for next update
  Omega <- matrix(0, nrow=nrow(Omega.tmp), ncol=ncol(Omega.tmp)) 
  
  # Will be referenced in update
  tmp <- Omega.tmp + t.step * (Sigma.tmp - S)

  for( a in unique(id) ){
    for(b in unique( id[id >= a] ) ) {
      a.ind <- which(id==a)
      b.ind <- which(id==b)
      
      new <- max(0, (1-t.step*lambda*
                       norm(tmp[a.ind, b.ind, drop=FALSE], type='F'))) *
            tmp[a.ind, b.ind]
      Omega[a.ind, b.ind] <- new
      if(max(abs(Omega[a.ind, b.ind] ))<1e-16) {Omega[a.ind, b.ind] <- 0}
      Omega[b.ind, a.ind] <- t(new)
      
    }
  }
  
  # Check if Omega is computationally PD
  keepPD <- isComputationallyPD(M=Omega) 
  if(det(Omega)>0){
    # Check that primal value increased
    new.primal <- sum(diag(S %*% Omega)) - log(det(Omega)) + 
      lambda*blockNorms(Omega, id, W) 
    old.primal <- sum(diag(S %*% Omega.tmp)) - log(det(Omega.tmp)) + 
      lambda*blockNorms(Omega.tmp, id, W) 
    keepPrimal <- TRUE 
  } else { keepPrimal <- FALSE } 
  
  if(keepPD & keepPrimal) {
    # Update covariance matrix
    Omega.new <- Omega
    s <- as.matrix(forceSymmetric(solve(Omega)))
    Sigma.new <- s
  
    if(det(Sigma.new)<0){
      # Undo it
      Omega.new <- Omega.tmp
      Sigma.new <- Sigma.tmp 
    }
    
  } else {
    Omega.new <- Omega.tmp
    Sigma.new <- Sigma.tmp
  }
    
  ans <- list(Omega=Omega.new, Sigma=Sigma.new)
  return(ans)
}


#' Generate lambdas to try
#' 
#' Generate penalty parameters spanning a range of 1 to 
#' \code{length(unique(id))} submatrix blocks.  This function is provided to 
#' establish a range of reasonable penalty parameters for optimization 
#' according to the algorithm of Kolar et al (2014).
#' 
#' @param S sample covariance
#' @param id vector of node identifiers
#' @param length.out number of lambda parameters to return
#' @return vector of \code{length.out} equally spaced lambdas
#' @export
#' 
#' @references
#' Kolar, M., Liu, H., and Xing, E. P. (2014). Graph estimation from 
#' multi-attribute data. The Journal of Machine Learning Research 15, 1713-1750.
#' 
#' @examples
#' Y <- matrix(rnorm(120), nrow=20, ncol=6)
#' S <- crossprod(Y)
#' id <- rep(1:3, each=2)
#' getLambdaRange(S, id, 5)
getLambdaRange <- function(S, id, length.out=10){
  fro <- blockNorms(M=S, id=id)
  min.fro <- min(fro[fro>0])
  diag(fro) <- 0
  fro <- as.numeric(fro)
  lambdas <- rev(seq(min.fro, max(fro), length.out=length.out+2))
  lambdas <- lambdas[2:(length.out+1)]
  return(lambdas)
}
