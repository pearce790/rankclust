#' Random data generation from Rank-Clustered BTL models
#'
#' This function randomly draws data from a Rank-Clustered Bradley-Terry-Luce (BTL) model with object-specific worth parameters, \code{omega}. Currently, only complete or partial rankings may be drawn, although those may be subsetted to simulate groupwise comparison data as a consequence of Luce's Axiom of Choice.
#'
#' @param I A numeric indicating the number of observations to draw.
#' @param J A numeric indicating the total number of objects being compared.
#' @param a_gamma A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.
#' @param b_gamma A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.
#' @param lambda A numeric for the Poisson hyperparameter on the number of non-empty clusters in the partition of worth parameters.
#' @param return_all A boolean indicator if K, partitions, nu, and omega should be return in addition to the generated rankings, Pi.
#'
#' @return A \code{I}x\code{R} matrix of rankings (if \code{return_all=FALSE}), or a list containing Pi, K, partitions, nu, and omega (if \code{return_all=TRUE}).
#'
#' @examples
#' set.seed(1)
#' rRCBTL(I=3,J=5,a_gamma=5,b_gamma=3,lambda=3,return_all=TRUE)
#' @export
rRCBTL <- function(I,J,a_gamma,b_gamma,lambda,return_all = FALSE){

  if((I %% 1)!=0 | I<1){stop("I must be a positive integer.")}
  if((J %% 1)!=0 | J<1){stop("J must be a positive integer.")}

  # sample number of unique clusters
  K <- sample(1:J,1,prob=dpois(1:J,lambda))

  # sample a partition given K
  size <- 0
  while(size < K){
    partition <- sample(1:K,J,replace=T)
    size <- length(unique(partition))
  }

  # sample nu_k, k=1,...,K
  nu <- rgamma(K,a_gamma,b_gamma)

  # determine omega using nu and partition
  omega <- nu[partition]

  # sample rankings
  Pi <- t(replicate(I,{sample(1:J,J,prob=omega)}))

  if(return_all){return(list(Pi=Pi,K=K,partition=partition,nu=nu,omega=omega))
  }else{return(Pi)}
}
