#' Partition sampling in Bayesian Rank-Clustered BTL models (internal use only)
#'
#' This function implements a reversible jump MCMC procedure for updating the parameter partition in Bayesian Rank-Clustered BTL models. For internal use only.
#'
#' @param Pi A matrix of preference orderings ("rankings"), such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumed that all unranked objects (if any) are less preferred than those which are ranked.
#' @param I A numeric indicating the number of rows in Pi
#' @param J A numeric indicating the total number of objects being compared.
#' @param nu A vector indicating current values for nu in the Gibbs sampler.
#' @param g A vector indicating current values for g in the Gibbs sampler.
#' @param K A vector indicating current values for K in the Gibbs sampler.
#' @param a_gamma A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.
#' @param b_gamma A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.
#' @param logprior_partition A J-vector indicating the log unnormalized probability for each possible partition, based on its number of non-empty clusters.
#' @param b_g The probability of "birth"ing a new partition cluster, if possible.
#' @param d_g The probability of "death"ing an existing partition cluster, if possible.
#' @param groupwise A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).
#'
#' @return A list containing updated values for g, nu, and K.
#'
#' @export
sample_partition <- function(Pi,I,J,nu,g,K,a_gamma,b_gamma,logprior_partition,b_g=0.5,d_g=0.5,groupwise=FALSE){
  S_g <- unlist(lapply(1:K,function(k){sum(g==k)}))
  if(rbinom(1,1,b_g)==1){
    ## Birth Mechanism
    #print("Birth!")
    # reject if all clusters have just one element
    if(all(S_g == 1)){
      logprob_accept <- -Inf
    }else{
      # pick a cluster to split
      which_split <- which(S_g>=2)
      if(length(which_split)==1){k <- which_split}else{k <- sample(which_split,1)}

      # decide how to split the objects
      new_clusts <- c(1,rbinom(S_g[k]-1,1,prob=0.5))
      while(sum(new_clusts)==0 | sum(new_clusts)==S_g[k]){
        new_clusts <- c(1,rbinom(S_g[k]-1,1,prob=0.5))
      }

      # get new values of nu and g
      u <- runif(1,min=0.5,max=1.5)
      nu_new <- nu
      nu_new[k] <- u*nu[k]
      nu_new[K+1] <- nu[k]/u
      g_new <- g
      g_new[which(g==k)[new_clusts==1]] <- K+1

      # calculation transition probability: check if nu_new are adjacent, else calculate A
      if(any((min(nu_new[c(k,K+1)]) <= nu[-k]) & (nu[-k] <= max(nu_new[c(k,K+1)])))){
        logprob_accept <- -Inf
      }else{
        logprob_accept <- dBTL(Pi=Pi,omega=nu_new[g_new],log=T,groupwise=groupwise)-dBTL(Pi=Pi,omega=nu[g],log=T,groupwise=groupwise)+
          sum(dgamma(nu_new[c(k,K+1)],a_gamma,b_gamma,log=T))-dgamma(nu[k],a_gamma,b_gamma,log=T)+
          logprior_partition[K+1]-logprior_partition[K]+
          log(d_g)+log(sum(S_g>=2))+log(2^(S_g[k]-1)-1)+log(2)+log(nu[k])-log(b_g)-log(K+1-1)-log(u)
      }
    }
    # accept/reject proposed g, nu, and K
    if(log(runif(1))<logprob_accept){
      # print("Accept!")
      # print(g_new)
      # print(nu_new)
      # print(nu_new[g_new])

      # sorting
      g <- unlist(lapply(1:J,function(j){which(sort(nu_new) == nu_new[g_new][j])}))
      nu <- sort(nu_new)
      K <- K+1
      if(any(nu_new[g_new] != nu[g])){print("Something wrong!")}

    }# else{print("Reject")}
    #print(nu[g])
    #print(g)
    #print(nu)
    #print(K)
  }else{
    ## Death Mechanism
    #print("Death!")
    # reject death if there is only one cluster
    if(K == 1){
      logprob_accept <- -Inf
    }else{
      # pick two adjacent clusters to merge
      if(K==2){which_merge <- c(1,2)}else{
        which_merge <- sample(1:(K-1),1)
        which_merge <- c(which_merge,which_merge+1)}

      # get new values of nu and g, and calculate u
      nu_new <- nu
      nu_new[which_merge[1]] <- sqrt(prod(nu[which_merge]))
      nu_new[which_merge[2]] <- NA
      g_new <- g
      g_new[g_new==which_merge[2]] <- which_merge[1]
      u <- nu_new[which_merge[1]]/nu[which_merge[1]]

      # calculate transition probability
      S_gnew <- unlist(lapply(1:max(g_new),function(k){sum(g_new==k)}))
      logprob_accept <- dBTL(Pi=Pi,omega=nu_new[g_new],log=T,groupwise=groupwise)-dBTL(Pi=Pi,omega=nu[g],log=T,groupwise=groupwise)+
        dgamma(nu_new[which_merge[1]],a_gamma,b_gamma,log=T)-sum(dgamma(nu[which_merge],a_gamma,b_gamma,log=T))+
        logprior_partition[K-1]-logprior_partition[K]+
        (log(b_g)+log(K-1)+log(u))-
        (log(d_g)+log(sum(S_gnew>=2))+log(2^(S_gnew[which_merge[1]]-1)-1)+log(2)+log(nu_new[which_merge[1]]))
    }
    # accept/reject proposed g, nu, and K
    if(log(runif(1))<logprob_accept){
      #print("Accept!")
      # print(g_new)
      # print(nu_new)
      #print(nu_new[g_new])

      #removing empty class and sorting
      g_new <- as.numeric(as.factor(g_new))
      nu_new <- c(na.exclude(nu_new))
      g <- unlist(lapply(1:J,function(j){which(sort(nu_new) == nu_new[g_new][j])}))
      nu <- sort(nu_new)
      K <- K-1
      if(any(nu_new[g_new] != nu[g])){print("Something wrong!")}

    }#else{print("Reject!")}
    #print(nu[g])
    #print(g)
    #print(nu)
    #print(K)
  }
  return(list(g=g,nu=nu,K=K))
}
