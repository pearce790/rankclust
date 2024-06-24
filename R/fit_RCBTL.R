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
