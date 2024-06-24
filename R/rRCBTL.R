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
