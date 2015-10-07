#' Multi-attribute subgraph estimation 
#' 
#' This estimates a single connected component from a multi-attribute graph. 
#' This function should only be called from within maLasso, as otherwise the
#' run will perform excessively large inversions.
#' 
#' @param S sample covariance matrix
#' @param n number of samples
#' @param id vector of node identifiers
#' @param lambda tuning parameter for penalty
#' @param W penalty weight matrix (same dimensions as S)
#' @param update how often to print updates
#' @param max.gap maximum allowable primal-dual gap
#' @param max.iter maximum number of iterations to complete (overriedes max.gap)
#' @param min.t minimum step size (overrides max.gap, max.iter)
#' @return list of precision, covariance, optimization status, lambda
#' 
#' @references
#' Kolar, M., Liu, H., and Xing, E. P. (2014). Graph estimation from 
#' multi-attribute data. The Journal of Machine Learning Research 15, 1713-1750.
#' 
#' @noRd
maLasso1 <- function(S, n, id, lambda, W, update=100, max.gap=0.5, max.iter=100, 
                     min.t=.Machine$double.eps){
  # Convenience calculations
  p <- ncol(S) # combined number of attributes
  nodes <- unique(id) # labels for nodes
  n.nodes <- length(unique(id)) # number of nodes
  k.nodes<- as.numeric(table(id)) # how many attributes for each
  
  # Initialize
  iter  <- 1
  t.step <- 2^-3
  if(ncol(S)>1){
    Sigma <- Omega <- diag(p)
  } else {
    Sigma <- S 
    Omega <- 1/S
  }

  primal <- primal.new <- getPrimal(Omega=Omega, S=S, W=W, id=id, lambda=lambda)
  dual <- getDual(Sigma)
  gap <- abs(primal-dual)
  soln <- NULL
  opt <- FALSE

  # Iterate
  while( iter < max.iter  & t.step >= min.t){
    # Propose a new precision matrix 
    Omega.new <- matrix(0, nrow=p, ncol=p)
    tmp <- Omega + t.step * (Sigma - S) 
    for(a in nodes){
      for(b in  nodes){
        a.ind <- which(id==a)
        b.ind <- which(id==b)
        Omega.new[a.ind,b.ind] <- max(0, 1-t.step*lambda/
                                norm(tmp[a.ind, b.ind, drop=FALSE], type='F')) *
          tmp[a.ind,b.ind, drop=FALSE]
        Omega.new[b.ind,a.ind] <- Omega.new[a.ind,b.ind]
      } 
    } 
    
    # Calculate a new primal value
    if( det(Omega) > 0 ){
      primal.new <- getPrimal(Omega=Omega.new, S=S, W=W, id=id, lambda=lambda)
    } 
    
    # If Omega is not computationally PD or if primal function has increased, 
    #   then decrase the step size and try the update again. Otherwise, update.
    if( !isComputationallyPD(Omega.new) | (primal.new >= primal) 
        | is.na(getDual(Sigma)) ){
      t.step <- t.step * 0.9
    } else {
      iter <- iter + 1 # only update iteration counter if update performed
      Omega <- Omega.new
      primal <- primal.new
      Sigma <- solve(Omega) # TODO: modify this later for larger problems
      Sigma <- (Sigma + t(Sigma)) / 2 # Force symmetry
      dual <- getDual(Sigma)
      gap <- abs(primal-dual)
      opt <- (gap < max.gap)
      
      # Update if requested 
      if(iter %% update == 0){
        print.str <- 'Iteration %d: primal=%.2E, dual=%.2E, gap=gap.now=%.2E, t=%.1E)'
        print(sprintf(print.str, iter, primal, dual, gap, t.step))
      }
    } 
    
    if(!is.na(getDual(Sigma)) &  abs(primal - getDual(Sigma) ) < max.gap ){
      opt <- TRUE 
      print.str <- 'CONVERGED at %d: gap =%.2E < %.2E, t=%.1E)'
      iter <- max.iter # get out of the loop
    } else if ( iter == max.iter | t.step < min.t | is.na(getDual(Sigma))) {
      print.str <- 'COMPLETED WITHOUT CONVERGENCE at %d: gap =%.2E > %.2E, t=%.1E)'
    }
  } # END: while( iter < max.iter  & t.step >= min.t)
  
  # Inform user whether the optimization was successful
  cat(sprintf(print.str, iter, gap, max.gap, t.step))
  
  # Assemble solution object
  soln <- list(Omega=Omega, Sigma=Sigma, opt=opt, lambda=lambda)
  return(soln)
}



#' Multi-attribute network estimation 
#' 
#' Estimates block precisions and covariances for a multi-attribute network
#' based on a Gaussian graphical model with zero mean vector.  This is according
#' to the method described by Kolar et al (2015).
#' 
#' @param S sample covariance matrix
#' @param n number of samples
#' @param id vector of node identifiers
#' @param lambda tuning parameter for penalty
#' @param W optional penalty weight matrix (same dimensions as S)
#' @param update how often to print updates
#' @param max.gap maximum allowable primal-dual gap
#' @param max.iter maximum number of iterations to complete (overriedes max.gap)
#' @param min.t minimum step size (overrides max.gap, max.iter)
#' @return list of precision, covariance, optimization status, lambda, and number components
#' 
#' @references
#' Kolar, M., Liu, H., and Xing, E. P. (2014). Graph estimation from 
#' multi-attribute data. The Journal of Machine Learning Research 15, 1713–1750.
#' 
#' @export
maLasso <- function(S, n, id, lambda, W=NULL, update=100, max.gap=0.5, 
                   max.iter=100, min.t=.Machine$double.eps){
  # Convenience calculations
  p <- ncol(S)
  nodes <- unique(id)
  n.nodes <- length(unique(id))
  
  # Check for actual weights
  if( is.null(W) ){ 
    # If NULL, then set as all entries = 1
    W <- matrix(1, nrow=n.nodes, ncol=n.nodes) 
  } else if ( any(dim(W)!=rep(n.nodes,2) )  ){ 
    # If supplied and incorrect dimensions, fail.
    stop('Error: dimensions of W incorrect')
  }
  
  # Find connected components
  fromat <- blockNorms(M=S, id=id)
  blockmat <- ifelse(abs(fromat) > lambda, 1, 0)
  block.graph <- graph.adjacency(blockmat, mode='undirected')
  block.comp <- clusters(block.graph)
  bl <- block.comp$membership
  cat(rbind(id, rep(bl, each=2)))
  
  # Build empty precision & covariance matrices
  Omega.tmp <- Sigma.tmp <- matrix(0, nrow=p, ncol=p)
  opt.tmp <- TRUE
  
  # Individual optimizations for each connected component
  for(b in sort(unique(bl))){
    cat(sprintf('Working on component %d of %d', b, length(unique((bl)))))
    which.nodes <- nodes[bl==b]
    which.cols <- id %in% which.nodes
    ans.bl <- maLasso1(S=S[which.cols, which.cols, drop=FALSE], n=n, 
                       id=id[which.cols],  lambda=lambda, 
                       W=W[which.nodes, which.nodes, drop=FALSE],
                       update=update, max.gap=max.gap, 
                       max.iter=max.iter, min.t=min.t)
    # Fill in values
    Omega.tmp[which.cols, which.cols] <- ans.bl$Omega
    Sigma.tmp[which.cols, which.cols] <- ans.bl$Sigma
    opt.tmp <- opt.tmp & ans.bl$opt
  }

  # Assemble solution object
  ans <- list(Omega=Omega.tmp, Sigma=Sigma.tmp, opt=opt.tmp, lambda=lambda, 
              n.comp=max(bl))
  return(ans)
}