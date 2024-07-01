#' Fit a Bayesian Rank-Clustered BTL model to ordinal comparison data (recommended for internal use only; use mcmc_RCBTL function instead)
#'
#' This function fits a Bayesian Rank-Clustered BTL model to ordinal comparison data (e.g., complete rankings, partial rankings, pairwise comparisons, or groupwise comparisons) such that the worth parameters are drawn from a PSSF prior (as defined in Pearce and Erosheva 2024).
#'
#' @param Pi A matrix of preference orderings ("rankings"), such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumed that all unranked objects (if any) are less preferred than those which are ranked.
#' @param J A numeric indicating the total number of objects being compared.
#' @param a_gamma A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.
#' @param b_gamma A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.
#' @param lambda A numeric for the Poisson hyperparameter on the number of non-empty clusters in the partition of worth parameters.
#' @param nu0 A numeric vector for the initialization of worth parameters, omega, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param num_iters A numeric indicating the total number of outer MCMC iterations (i.e., the number of times the partition is updated in the Gibbs sampler).
#' @param nu_reps A numeric indicating the number of times each worth parameter is drawn per update of the parameter partition. There will be a total of \code{num_iters}x\code{nu_reps} samples from the posterior.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#'
#' @return A list of 4 elements: \code{omega}, a (\code{num_iters}x\code{nu_reps})x\code{J} matrix of approximate posterior draws of the object-specific worth parameters, omega; \code{nu} a (\code{num_iters}x\code{nu_reps})x\code{J} matrix of the unique parameter values corresponding to the jth partition cluster in posterior draw i, \code{g} a a (\code{num_iters}x\code{nu_reps})x\code{J} matrix indicating the cluster membership of object j in posterior draw i, and \code{K} a vector of the number of non-empty partition clusters in each posterior draw.
#'
#' @examples
#' Pi <- matrix(data=c(1,2,3,NA,NA,1,2,3,4,5),byrow=TRUE,nrow=2)
#' fit_RCBTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,lambda=2,num_iters=5,nu_reps=2)
#' @export
fit_RCBTL <- function(Pi,J,a_gamma,b_gamma,lambda,nu0=NULL,num_iters=100,nu_reps=3,groupwise=FALSE){

  ### Step 0 ###

  # calculate constants and fixed values
  Pi <- as.matrix(Pi)
  I <- nrow(Pi)
  logprior_partition <- dpois(1:J,lambda=lambda,log=T)
  constants <- obtain_constants(Pi=Pi,I=I,J=J,groupwise=groupwise)

  # setup storage for gibbs samples
  nu_samples <- matrix(NA,nrow=num_iters*nu_reps,ncol=J)
  g_samples <- matrix(NA,nrow=0,ncol=J)
  K_samples <- c()
  curr <- 1

  ### Step 1 ###
  if(is.null(nu0)){
    K <- sample(1:J,1,prob=exp(logprior_partition))
    nu <- rgamma(K,a_gamma,b_gamma)
    if((J-K) > 0){g <- c(1:K,sample(1:K,J-K,replace=T))}else{g <- 1:K}
  }else{
    g <- as.numeric(factor(nu0))
    K <- length(unique(g))
    nu <- as.numeric(levels(factor(nu0)))
  }

  ### Step 2 ###
  for(iter in 1:num_iters){
    ### Step 2(a): Sample new partition, g (and update nu and K accordingly) ###
    g_curr <- sample_partition(nu=nu,g=g,K=K,Pi=Pi,I=I,J=J,a_gamma=a_gamma,b_gamma=b_gamma,
                               logprior_partition = logprior_partition,groupwise=groupwise)
    g <- g_curr$g
    nu <- g_curr$nu
    K <- g_curr$K

    ### Step 2(b): Sample new worth parameters, nu (and store associated nu and K) ###
    nu_curr <- sample_nu(n=nu_reps,Pi=Pi,I=I,J=J,nu=nu,p=g,K=K,a_gamma=a_gamma,b_gamma=b_gamma,
                         c_k0=constants$c_k0,delta_irk0=constants$delta_irk0,groupwise=groupwise)
    nu_samples[curr:(curr+nu_reps-1),1:K] <- rbind(nu_curr)
    g_samples <- rbind(g_samples,matrix(rep(g,nu_reps),ncol=J,byrow=T))
    K_samples <- c(K_samples,rep(K,nu_reps))

    # housekeeping: keep track of final "nu" and update counter
    nu <- nu_curr[nu_reps,]
    curr <- curr+nu_reps
  }

  # calculate omega for each iteration
  omega_samples <- matrix(unlist(lapply(1:nrow(nu_samples),function(iter){
    nu_samples[iter,g_samples[iter,]]
  })),nrow=nrow(nu_samples),byrow=T)

  # return results
  return(list(omega=omega_samples,nu=nu_samples,g=g_samples,K=K_samples))
}
