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
