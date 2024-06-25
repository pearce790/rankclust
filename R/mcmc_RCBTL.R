#' Fit Bayesian Rank-Clustered BTL model to ordinal comparison data using multiple MCMC chains
#'
#' This function fits a Bayesian Rank-Clustered BTL model to ordinal comparison data (e.g., complete rankings, partial rankings, pairwise comparisons, or groupwise comparisons) such that the worth parameters are drawn from a PSSF prior (as defined in Pearce and Erosheva 2024). The function has input parameters to permit drawing multiple MCMC chains, as well as thinning and burn-in.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @import utils
#'
#' @param Pi A matrix of rankings, such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumd that all unranked objects (if any) are less preferred than those which are ranked.
#' @param J A numeric indicating the total number of objects being compared.
#' @param a_gamma A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.
#' @param b_gamma A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.
#' @param lambda A numeric for the Poisson hyperparameter on the number of non-empty clusters in the partition of worth parameters.
#' @param nu0 A numeric vector for the initialization of worth parameters, omega, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param num_iters A numeric indicating the total number of outer MCMC iterations (i.e., the number of times the partition is updated in the Gibbs sampler).
#' @param nu_reps A numeric indicating the number of times each worth parameter is drawn per update of the parameter partition. There will be a total of \code{num_iters}x\code{nu_reps} samples from the posterior.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#' @param chains A numeric indicating the total number of independent MCMC chains to be run.
#' @param burn_prop A numeric between 0 and 1 indicating the proportion of MCMC samples in each chain to be removed as burn-in.
#' @param thin A numeric indicating that only every \code{thin}-th sample should be retained, to save computational memory.
#' @param seed A numeric indicating the random seed that should be set before running the first MCMC chain.
#' @param normalize_omega A boolean indicating if each posterior draw of omega should be normalized post-hoc to sum to 1; removes a standard identifiability concern of BTL models.
#'
#' @return A (\code{chains}x\code{num_iters}/\code{thin})x(3J+3) matrix of posterior draws, one row per posterior sample of omega, nu, and g, with additional columns indicating the MCMC chain index, iteration index, and number of non-empty partition clusters K of each posterior sample.
#'
#' @examples
#' Pi <- matrix(data=c(1,2,3,NA,NA,1,2,3,4,5),byrow=TRUE,nrow=2)
#' mcmc_RCBTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,lambda=2,num_iters=6,chains=2,seed=1)
#' @export
mcmc_RCBTL <- function(Pi,J,a_gamma,b_gamma,lambda,nu0=NULL,num_iters=1000,nu_reps=2,groupwise=FALSE,
                       chains=4,burn_prop=0.5,thin=1,seed=NULL,normalize_omega=TRUE){
  if(!is.null(seed)){set.seed(seed)}
  counter <- 1
  mcmc <- replicate(n=chains,{
    print(paste0("Estimating chain ",counter," of ",chains,"."))
    counter <<- counter+1
    res <- fit_RCBTL(Pi=Pi,J=J,a_gamma=a_gamma,b_gamma=b_gamma,lambda=lambda,nu0=nu0,
                     num_iters=num_iters,nu_reps=nu_reps,groupwise=groupwise)
    nreps <- nrow(res$omega)
    keep_reps <- seq(ceiling(burn_prop*nreps)+1,nreps,by=thin)

    if(normalize_omega){
      omega <- t(apply(res$omega[keep_reps,,drop=FALSE],1,function(x){x/sum(x)}))
    }else{omega <- res$omega[keep_reps,,drop=FALSE]}
    tmp <- as.data.frame(cbind(omega,res$nu[keep_reps,,drop=FALSE],
                               res$g[keep_reps,,drop=FALSE],res$K[keep_reps],keep_reps))
    names(tmp) <- c(paste0("omega",1:J),paste0("nu",1:ncol(res$nu)),paste0("G",1:J),"K","iteration")
    return(tmp)
  },simplify = FALSE)
  mcmc <- do.call(rbind,mcmc)
  mcmc$chain <- factor(rep(1:chains,each=nrow(mcmc)/chains))
  mcmc <- mcmc %>% dplyr::select(chain,iteration,K,dplyr::everything())
  rm(counter)

  return(mcmc)
}
