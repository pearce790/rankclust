#' Calculate the density of ordinal comparison data based on BTL models
#'
#' This function calculates the density of ordinal comparison data (e.g., complete rankings, partial rankings, pairwise comparisons, or groupwise comparisons) based on an appropriate Bradley-Terry-Luce (BTL) model with object-specific worth parameters, \code{omega}.
#'
#' @import stats
#'
#' @param Pi A matrix of preference orderings ("rankings"), such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumed that all unranked objects (if any) are less preferred than those which are ranked.
#' @param omega A vector of non-negative object worth parameters.
#' @param log A boolean to indicate whether the log density should be returned. Default to \code{FALSE}.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#'
#' @return A numeric joint density of the observed rankings.
#'
#' @examples
#' Pi <- matrix(data=c(1,2,3,NA,NA,1,2,3,4,5),byrow=TRUE,nrow=2)
#' dBTL(Pi=Pi,omega=c(10,8,6,4,2),log=TRUE,groupwise=FALSE)
#' dBTL(Pi=Pi,omega=c(10,8,6,4,2),log=TRUE,groupwise=TRUE)
#' @export
dBTL <- function(Pi,omega,log=FALSE,groupwise=FALSE){

  # ensure proper data format
  if(!is.matrix(Pi)){stop("Pi must be a matrix with one row per observation")}
  omega <- c(omega)
  if(!is.logical(log)){stop("log must be logical value (TRUE/FALSE)")}

  # calculate J,theta
  J <- length(omega)
  theta <- log(omega)

  # calculate probability terms
  if(groupwise){
    terms_log <- t(apply(Pi,1,function(pi){
      Ri <- length(na.exclude(pi))
      if(Ri == (J-1)){stop("Pi invalid: Cannot have Ri = J-1")}
      terms <- unlist(lapply(1:Ri,function(r){theta[pi[r]]-matrixStats::logSumExp(theta[setdiff(pi[1:Ri],pi[0:(r-1)])])}))
      c(terms,rep(NA,J-Ri))
    }))
  }else{
    terms_log <- t(apply(Pi,1,function(pi){
      Ri <- length(na.exclude(pi))
      if(Ri == (J-1)){stop("Pi invalid: Cannot have Ri = J-1")}
      terms <- unlist(lapply(1:Ri,function(r){theta[pi[r]]-matrixStats::logSumExp(theta[setdiff(1:J,pi[0:(r-1)])])}))
      c(terms,rep(NA,J-Ri))
    }))
  }

  # return result
  if(log){return(sum(terms_log,na.rm=T))
  }else{return(exp(sum(terms_log,na.rm=T)))}
}
