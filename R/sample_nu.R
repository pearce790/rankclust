sample_nu <- function(n=1,Pi,I,J,nu,p,K,a_gamma,b_gamma,c_k0=NULL,delta_irk0=NULL,groupwise=FALSE){

  # calculate values c_k and delta_irk
  c_k <- rep(NA,K)
  delta_irk <- array(NA,c(I,J,K))
  Ri <- apply(Pi,1,function(pi){length(na.exclude(pi))})
  if(is.null(c_k0) | is.null(delta_irk0)){
    # calculate values from scratch, if helper constants are unavailable (slower!)
    flat_Pi <- na.exclude(c(Pi))
    for(k in 1:K){
      which_k <- which(p==k)
      c_k[k] <- sum(flat_Pi %in% which_k)
      for(i in 1:I){
        pi_filled <- Pi[i,1:Ri[i]]
        if(groupwise==FALSE){pi_filled <- c(pi_filled,setdiff(1:J,pi_filled))}
        vals <- pi_filled %in% which_k
        delta_irk[i,1:Ri[i],k] <- unlist(lapply(1:Ri[i],function(r){sum(vals[r:J],na.rm=T)}))
      }
    }
    c_k_slow <- c_k
    delta_irk_slow <- delta_irk
  }else{
    # calculate values for current partition p and K if helper constants have been pre-calculated (faster!)
    for(k in 1:K){
      which_k <- which(p==k)
      c_k[k] <- sum(c_k0[which_k])
      delta_irk[,,k] <- apply(delta_irk0[,,which_k],c(1,2),sum)
    }
  }

  # update nu
  nu_new <- matrix(data=NA,nrow=n,ncol=K)
  for(iter in 1:n){
    # sample auxiliary variables, Y_ir
    Y <- matrix(NA,nrow=I,ncol=J)
    for(i in 1:I){
      if(groupwise==FALSE){v_curr <- sum(nu[p])
      }else{v_curr <- sum(nu[p][Pi[i,]],na.rm=T)}
      for(r in 1:Ri[i]){
        Y[i,r] <- rexp(1,rate=v_curr)
        v_curr <- v_curr - nu[p][Pi[i,r]]
      }
    }

    # sample updated nu
    nu_new[iter,] <- nu <- unlist(lapply(1:K,function(k){rgamma(1,a_gamma+c_k[k],b_gamma+sum(Y * delta_irk[,,k],na.rm=T))}))
  }

  return(nu_new)
}
