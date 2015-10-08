#' Wrapper function for multi-attribute network estimation via lasso
#' 
#' Performs estimation of a zero-mean multi-attribute Gaussian graphical model 
#' as described by Kolar et al (2014) when \code{mode=1}.
#' 
#' If \code{mode=2}, the \code{id} argument is essentially disregarded, and all
#' rows/columns in \code{S} are treated as if they represent individual nodes. 
#' If \code{mode=3}, the optimization of \code{mode=2} is performed separately
#' for each attribute type, and a matrix with 0 on all cross-attribute entries
#' is returned.
#' 
#' @param S sample covariance matrix
#' @param n number of samples
#' @param id vector assigning variables to nodes
#' @param lambda penalty tuning parameter
#' @param W optional weight matrix 
#' @param mode estimation mode. Can be 1 (multiattribute; default), 
#' 2 (unstructured together), or 3 (separately; only valid for 2-attribute data)
#' @param update how often to print updates in optimization
#' @param max.gap maximum allowable primal/dual gap
#' @param max.iter maximum number of iterations to optimize (overrides max.gap)
#' @param min.t minimum step size (overrides max.gap, max.iter)
#' @return list of precision, covariance, optimization status, lambda, 
#' and number components
#' 
#' @references
#' Kolar, M., Liu, H., and Xing, E. P. (2014). Graph estimation from 
#' multi-attribute data. The Journal of Machine Learning Research 15, 1713-1750.
#' 
#' @examples
#' library(mvtnorm)
#' id <- rep(1:3, each=2)
#' Omega <- matrix(0, nrow=6, ncol=6)
#' Omega[1:4,1:4] <- 1
#' diag(Omega) <- 6
#' Sigma <- solve(Omega)
#' n <- 1000
#' Y <- rmvnorm(n, sigma=Sigma)
#' S <- crossprod(Y)
#' 
#' lambdas <- getLambdaRange(S, id, 10)
#' result <- multiAttEstimate(S, n, id, lambdas[1])
#' 
#' @export
multiAttEstimate <- function(S, n, id, lambda, W=NULL, mode=1, update=100, 
                     max.gap=0.5, max.iter=100, min.t=.Machine$double.eps){
  # Mode can be 1 (multiattribute), 2 (unstructured together), 3 (separately)
  if(mode==1) {
    # Run multiattribute
    cat('Running in multi-attribute mode\n')
    ans <- multiAttLasso(S=S, n=n, id=id, lambda=lambda, W=W, update=update, 
                    max.gap=max.gap, max.iter=max.iter, min.t=min.t)
    
  } else if (mode==2) {
    # Run unstructured together
    cat('Running in unstructured mode\n')
    id2 <- 1:ncol(S)
    ans <- multiAttLasso(S=S, n=n, id=id2, lambda=lambda, W=W, update=update, 
                    max.gap=max.gap, max.iter=max.iter, min.t=min.t)
    
  } else if (mode==3) {
    # Run separately 
    cat('Running in separated mode\n')
    id3a <- which(!duplicated(id))
    id3b <- which(duplicated(id))
    cat('First optimization...\n')
    ans3a <- multiAttLasso(S=S[id3a,id3a,drop=FALSE], n=n, id=1:length(id3a), 
                     lambda=lambda, W=W,  update=update, max.gap=max.gap, 
                     max.iter=max.iter, min.t=min.t)
    
    cat('Second optimization...\n')
    ans3b <- multiAttLasso(S=S[id3b,id3b,drop=FALSE],n=n, id=1:length(id3b), 
                     lambda=lambda, W=W,  update=update, max.gap=max.gap, 
                     max.iter=max.iter, min.t=min.t)
    
    Omega.tmp <- Sigma.tmp <- matrix(0, nrow=ncol(S), ncol=ncol(S))
    Omega.tmp[id3a, id3a] <- ans3a$Omega
    Omega.tmp[id3b, id3b] <- ans3b$Omega
    
    Sigma.tmp[id3a, id3a] <- ans3a$Sigma
    Sigma.tmp[id3b, id3b] <- ans3b$Sigma
    
    opt.tmp <- ans3a$opt & ans3b$opt
    
    ans <- list(Omega=Omega.tmp, Sigma=Sigma.tmp, opt=opt.tmp, lambda=lambda, 
                n.comp=ans3a$n.comp+ans3b$n.comp)
  } else {
    stop('Error: mode must be 1, 2, or 3') 
  }
  
  # Return whatever solution you've generated
  return(ans)
}


#' Wrapper function for multi-attribute network estimation plus selection via lasso
#' 
#' Performs estimation of a zero-mean multi-attribute Gaussian graphical model 
#' as described by Kolar et al (2014) when \code{mode=1} with a selection 
#' procedure to determine optimal \code{lambda}.  Selection may be performed
#' according to extended BIC, BIC, AIC, or Hamming distance to true graph.
#' The true graph \code{Theta} is only required if Hamming distance is being
#' used for the selection procedure. 
#' 
#' If \code{mode=2}, the \code{id} argument is essentially disregarded, and all
#' rows/columns in \code{S} are treated as if they represent individual nodes. 
#' If \code{mode=3}, the optimization of \code{mode=2} is performed separately
#' for each attribute type, and a matrix with 0 on all cross-attribute entries
#' is returned.
#' 
#' @param S sample covariance matrix
#' @param n number of samples
#' @param id vector assigning variables to nodes
#' @param lambda.range vector of lambdas to try
#' @param W optional weight matrix
#' @param mode estimation mode. Can be 1 (multiattribute; default), 
#' 2 (unstructured together), or 3 (separately; only valid for 2-attribute data)
#' @param update how often to print updates in optimization
#' @param max.gap maximum allowable primal/dual gap
#' @param max.iter maximum number of iterations to optimize (overrides max.gap)
#' @param min.t minimum step size (overrides max.gap, max.iter
#' @param method how to select the best model (\code{'EBIC'}, \code{'BIC'}, 
#' \code{'AIC'}, or \code{'hamming'})
#' @param Theta.true true underlying graph (\code{'hamming'} method only)
#' @param plot boolean, whether or not to make diagnostic plot
#' @param gamma gamma parameter for EBIC method (default 0.5)
#' @return list of precision, covariance, optimization status, lambda, and 
#' number components
#' 
#' @references
#' Kolar, M., Liu, H., and Xing, E. P. (2014). Graph estimation from 
#' multi-attribute data. The Journal of Machine Learning Research 15, 1713-1750.
#' 
#' Chen, J. and Chen, Z. (2008). Extended Bayesian information criteria for 
#' model selection with large model spaces. Biometrika 95, 759-771.
#' 
#' @examples
#' library(mvtnorm)
#' id <- rep(1:3, each=2)
#' Omega <- matrix(0, nrow=6, ncol=6)
#' Omega[1:4,1:4] <- -1
#' diag(Omega) <- 6
#' Sigma <- solve(Omega)
#' n <- 1000
#' Y <- rmvnorm(n, sigma=Sigma)
#' S <- crossprod(Y)
#' 
#' lambdas <- getLambdaRange(S, id, 3)
#' result <- multiAttSelect(S, n, id, lambdas)
#' 
#' 
#' @export
multiAttSelect <- function(S,n, id, lambda.range, W=NULL, mode=1, update=100, 
                           max.gap=0.5, max.iter=100, min.t=.Machine$double.eps,
                           method='BIC', Theta.true=NULL, plot=NULL, 
                           gamma=0.5){  
  # Printing
  if(method=='hamming' & is.null(Theta.true)){
    cat('Error: must specify Theta.true if selecting by Hamming distance.  
            Will run in BIC mode instead.\n')
    method <- 'BIC'
  }
  if(method=='BIC'){
    cat(sprintf('----- Finding lambda by min BIC: lambda=(%s)\n', 
                  paste(lambda.range, collapse=', ')))
  } else if(method=='AIC'){
    cat(sprintf('----- Finding lambda by min AIC: lambda=(%s)\n', 
                  paste(lambda.range, collapse=', ')))
  } else if(method=='EBIC'){
    cat(sprintf('----- Finding lambda by min EBIC (gamma=%.2f): lambda=(%s)\n', 
                  gamma, paste(lambda.range, collapse=', ')))
  }
  else if (method=='hamming'){
    cat(sprintf('----- Finding lambda by min Hamming distance: lambda=(%s)\n', 
                  paste(lambda.range, collapse=', ')))
  }

  # Setup
  p <- ncol(S)
  tuning <- data.frame(lambda=lambda.range, n.comp=NA, BIC=NA, EBIC=NA, AIC=NA, 
                       opt=NA)
  if(!is.null(Theta.true)){
    tuning$hamming <- NA
  }

  # For each value of lambda, run.  Save ans & BIC if first or best
  for(lambda in lambda.range){
    cat(sprintf('+++++ Trying lambda=%.4f', lambda))
    # Run optimization
    ans.now <- multiAttEstimate(S=S, n=n, id=id, lambda=lambda, W=W, mode=mode, 
                        update=update, max.gap=max.gap, max.iter=max.iter, 
                        min.t=min.t)
    if(mode==1){
      id.use <- id
    } else {
      id.use <- 1:p 
    }
      
    BIC.now <- getBIC(S=S, n=n, Omega=ans.now$Omega,id=id.use)
    EBIC.now <- getEBIC(S=S, n=n, Omega=ans.now$Omega,id=id.use, gamma=gamma)
    AIC.now <- getAIC(S=S, n=n, Omega=ans.now$Omega,id=id.use)
    opt.now <- ans.now$opt
    if(!is.null(Theta.true)){
      hamming.now <- getHammingDist(Omega=ans.now$Omega, Theta=Theta.true, 
                                    id=id)
      tuning[tuning$lambda==lambda, 'hamming'] <- hamming.now
    }
    
    # Depending on method setting, set performance parameter
    if(method=='BIC'){
      measure.now <- BIC.now
    } else if (method=='AIC'){
      measure.now <- AIC.now
    } else if (method=='hamming'){
      measure.now <- hamming.now
    } else if (method=='EBIC'){
      measure.now <- EBIC.now
    }
    
    
    # Fill in tuning data frame
    tuning[tuning$lambda==lambda, method] <- measure.now
    tuning[tuning$lambda==lambda, 'AIC'] <- AIC.now
    tuning[tuning$lambda==lambda, 'BIC'] <- BIC.now
    tuning[tuning$lambda==lambda, 'EBIC'] <- EBIC.now
    tuning[tuning$lambda==lambda,'opt'] <- opt.now
    tuning[tuning$lambda==lambda, 'n.comp'] <- ans.now$n.comp
    
    # If it's the smallest BIC, then keep the answer
    # If not optimal and ans$opt
    if( lambda==lambda.range[1] ){
      # If this is the first time through loop
      ans <- ans.now
      lambda.best <- lambda
    } else if(measure.now == min(tuning[,method], na.rm=TRUE)){
      if( BIC.now==min(tuning$BIC[tuning[,method]==measure.now], na.rm=TRUE)){
        # If there's a tie (method 'hamming'), break via BIC
        ans <- ans.now
        lambda.best <- lambda
      } 
    }
  }
  tuning$Best <- ifelse(tuning$lambda == lambda.best, 'BEST', ' ')
  cat(sprintf('----- BEST %s=%.2f at lambda==%.2f\n', method, 
                min(tuning[,method]), lambda.best))
  print(tuning)

  cols <- c('lambda','n.comp','AIC','BIC','EBIC','hamming')
  cols <- cols[cols %in% names(tuning)]
  if(!is.null(plot)){
    png(sprintf('tuning_%s_mode%d.png', plot, mode))
    pairs(tuning[,cols],
         col=ifelse(tuning$lambda==lambda.best, 'red','black'), 
          main=sprintf('%s, mode %d', plot, mode), pch=16)
    dev.off()
  }

  return(ans)
}