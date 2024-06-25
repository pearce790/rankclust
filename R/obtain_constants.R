#' Calculation of constants for Gibbs sampling of worth parameters via data augmentation (internal use only)
#'
#' This function calculates constants c_k0 and delta_irk0 based on the observed rankings \code{Pi}, as defined in Pearce and Erosheva (2024). For internal use only.
#'
#' @import stats
#'
#' @param Pi A matrix of rankings, such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumd that all unranked objects (if any) are less preferred than those which are ranked.
#' @param I A numeric indicating the number of rows in Pi
#' @param J A numeric indicating the total number of objects being compared.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#'
#' @return A list containing calculated constants c_k0 and delta_irk0.
#'
#' @export
obtain_constants <- function(Pi,I,J,groupwise=FALSE){
  c_k0 <- rep(NA,J)
  delta_irk0 <- array(NA,c(I,J,J))
  flat_Pi <- na.exclude(c(Pi))
  Ri <- apply(Pi,1,function(pi){length(na.exclude(pi))})
  for(k in 1:J){
    which_k <- which(1:J==k)
    c_k0[k] <- sum(flat_Pi %in% which_k)
    for(i in 1:I){
      pi_filled <- Pi[i,1:Ri[i]]
      if(groupwise==FALSE){pi_filled <- c(pi_filled,setdiff(1:J,pi_filled))}
      vals <- pi_filled %in% which_k
      delta_irk0[i,1:Ri[i],k] <- unlist(lapply(1:Ri[i],function(r){sum(vals[r:J],na.rm=T)}))
    }
  }
  return(list(c_k0=c_k0,delta_irk0=delta_irk0))
}
