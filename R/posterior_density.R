#' Calculate proportional posterior density for mcmc samples in a Rank-Clustered BTL model.
#'
#' This function calculates the (proportional) posterior density for MCMC samples from a Rank-Clustered BTL model.
#'
#' @import stats
#'
#' @param Pi A matrix of preference orderings ("rankings"), such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumed that all unranked objects (if any) are less preferred than those which are ranked.
#' @param mcmc The matrix output of the \code{mcmc_BTL} or \code{mcmc_RCBTL} functions.
#' @param a_gamma A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.
#' @param b_gamma A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.
#' @param lambda A numeric for the Poisson hyperparameter on the number of non-empty clusters in the partition of worth parameters.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#' @param log A boolean to indicate whether the log posterior densities should be returned. Default to \code{FALSE}.
#'
#' @return A ggplot(s) of side-by-side violin plots for each omega_j.
#'
#' @examples
#' Pi <- rBTL(I=30,omega=c(4^2,4^2,4^1,4^1,4^0))
#' mcmc <- mcmc_RCBTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,lambda=1,num_iters=200,chains=2,seed=1)
#' posterior_density(Pi = Pi, mcmc = mcmc, a_gamma=5, b_gamma = 3, lambda=2, log=TRUE)
#' @export
posterior_density <- function(Pi,mcmc,a_gamma,b_gamma,lambda,groupwise=FALSE,log=FALSE){
  log_K <- dpois(mcmc$K,lambda=lambda,log=T)
  log_nu <- apply(mcmc[,grep("nu",names(mcmc))],1,function(nu_values){sum(dgamma(nu_values,a_gamma,b_gamma,log=T),na.rm=T)})
  log_Pi <- apply(mcmc[,grep("omega",names(mcmc))],1,function(omega){dBTL(Pi,omega,log=T,groupwise=groupwise)})
  log_density <- log_K + log_nu + log_Pi

  if(log){return(log_density)
  }else{return(exp(log_density))}
}
