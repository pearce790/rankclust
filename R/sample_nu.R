#' Sampling of unique parameter worth values, nu, in Bayesian Rank-Clustered BTL models (internal use only)
#'
#' This function implements a Gibbs sampler via data augmentation for updating the unique worth parameters, nu, in Bayesian Rank-Clustered BTL models. For internal use only.
#'
#' @param n A numeric indicating the number of posterior samples to draw.
#' @param Pi A matrix of rankings, such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumd that all unranked objects (if any) are less preferred than those which are ranked.
#' @param I A numeric indicating the number of rows in Pi
#' @param J A numeric indicating the total number of objects being compared.
#' @param nu A vector indicating current values for nu in the Gibbs sampler.
#' @param p A vector indicating current values for g in the Gibbs sampler.
#' @param K A vector indicating current values for K in the Gibbs sampler.
#' @param a_gamma A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.
#' @param b_gamma A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.
#' @param c_k0 A vector of constants, as calculated by the obtain_constants function.
#' @param delta_irk0 A matrix of constants, as calculated by the obtain_constants function.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#'
#' @return A matrix containing updated values for nu for n posterior samples.
#'
#' @export
sample_nu <- function(n=1,Pi,I,J,nu,p,K,a_gamma,b_gamma,c_k0=NULL,delta_irk0=NULL,groupwise=FALSE){

  # calculate values c_k and delta_irk
  c_k <- rep(NA,K)
  delta_irk <- array(NA,c(I,J,K))
  Ri <- apply(Pi,1,function(pi){length(na.exclude(pi))})
  if(is.null(c_k0) | is.null(delta_irk0)){
    # calculate values from scratch, if helper constants are unavailable (slower!)
    flat_Pi <- na.exclude(c(Pi))
    for(k in 1:K){
      which_k <- which(p==k)
      c_k[k] <- sum(flat_Pi %in% which_k)
      for(i in 1:I){
        pi_filled <- Pi[i,1:Ri[i]]
        if(groupwise==FALSE){pi_filled <- c(pi_filled,setdiff(1:J,pi_filled))}
        vals <- pi_filled %in% which_k
        delta_irk[i,1:Ri[i],k] <- unlist(lapply(1:Ri[i],function(r){sum(vals[r:J],na.rm=T)}))
      }
    }
    c_k_slow <- c_k
    delta_irk_slow <- delta_irk
  }else{
    # calculate values for current partition p and K if helper constants have been pre-calculated (faster!)
    for(k in 1:K){
      which_k <- which(p==k)
      c_k[k] <- sum(c_k0[which_k])
      delta_irk[,,k] <- apply(delta_irk0[,,which_k],c(1,2),sum)
    }
  }

  # update nu
  nu_new <- matrix(data=NA,nrow=n,ncol=K)
  for(iter in 1:n){
    # sample auxiliary variables, Y_ir
    Y <- matrix(NA,nrow=I,ncol=J)
    for(i in 1:I){
      if(groupwise==FALSE){v_curr <- sum(nu[p])
      }else{v_curr <- sum(nu[p][Pi[i,]],na.rm=T)}
      for(r in 1:Ri[i]){
        Y[i,r] <- rexp(1,rate=v_curr)
        v_curr <- v_curr - nu[p][Pi[i,r]]
      }
    }

    # sample updated nu
    nu_new[iter,] <- nu <- unlist(lapply(1:K,function(k){rgamma(1,a_gamma+c_k[k],b_gamma+sum(Y * delta_irk[,,k],na.rm=T))}))
  }

  return(nu_new)
}
