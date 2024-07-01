#' Fit a Bayesian BTL model to ordinal comparison data (recommended for internal use only; use mcmc_BTL function instead)
#'
#' This function fits a standard Bayesian BTL model to ordinal comparison data (e.g., complete rankings, partial rankings, pairwise comparisons, or groupwise comparisons) such that each worth parameter receives an independent Gamma prior.
#'
#' @param Pi A matrix of preference orderings ("rankings"), such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumed that all unranked objects (if any) are less preferred than those which are ranked.
#' @param J A numeric indicating the total number of objects being compared.
#' @param a_gamma A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.
#' @param b_gamma A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.
#' @param nu0 A numeric vector for the initialization of worth parameters, omega, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param num_iters A numeric indicating the total number of MCMC iterations.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#'
#' @return A list containing a single element, \code{omega}, a \code{num_iters}x\code{J} matrix of approximate posterior draws of the object-specific worth parameters, omega.
#'
#' @examples
#' Pi <- matrix(data=c(1,2,3,NA,NA,1,2,3,4,5),byrow=TRUE,nrow=2)
#' fit_BTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,num_iters=10)
#' @export
fit_BTL <- function(Pi,J,a_gamma,b_gamma,nu0=NULL,num_iters=100,groupwise=FALSE){

  ### Step 0 ###

  # calculate constants and fixed values
  Pi <- as.matrix(Pi)
  I <- nrow(Pi)
  constants <- obtain_constants(Pi=Pi,I=I,J=J,groupwise=groupwise)

  # setup storage for gibbs samples
  nu_samples <- matrix(NA,nrow=num_iters,ncol=J)
  curr <- 1

  ### Step 1 ###
  if(is.null(nu0)){nu <- rgamma(J,a_gamma,b_gamma)
  }else{nu <- nu0}
  g <- 1:J

  ### Step 2 ###
  for(iter in 1){
    ### Step 2(a) ###
    # not applicable for regular BTL model

    ### Step 2(b): Sample new worth parameters, nu (and store associated nu and K) ###
    nu_curr <- sample_nu(n=num_iters,Pi=Pi,I=I,J=J,nu=nu,p=g,K=J,a_gamma=a_gamma,b_gamma=b_gamma,
                         c_k0=constants$c_k0,delta_irk0=constants$delta_irk0,groupwise=groupwise)
    nu_samples[curr:(curr+num_iters-1),1:J] <- rbind(nu_curr)
  }



  # return results
  return(list(omega=nu_samples))
}
