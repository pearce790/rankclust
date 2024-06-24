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
