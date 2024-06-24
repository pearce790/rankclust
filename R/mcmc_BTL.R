mcmc_BTL <- function(Pi,J,a_gamma,b_gamma,nu0=NULL,num_iters=1000,groupwise=FALSE,
                     chains=4,burn_prop=0.5,thin=1,seed=NULL,normalize_omega=TRUE){

  if(!is.null(seed)){set.seed(seed)}
  counter <- 1
  mcmc <- replicate(n=chains,{
    print(paste0("Estimating chain ",counter," of ",chains,"."))
    counter <<- counter+1
    res <- fit_BTL(Pi=Pi,J=J,a_gamma=a_gamma,b_gamma=b_gamma,nu0=nu0,num_iters=num_iters,groupwise=groupwise)
    nreps <- nrow(res$omega)
    keep_reps <- seq(ceiling(burn_prop*nreps)+1,nreps,by=thin)

    if(normalize_omega){
      omega <- t(apply(res$omega[keep_reps,,drop=FALSE],1,function(x){x/sum(x)}))
    }else{omega <- res$omega[keep_reps,,drop=FALSE]}
    tmp <- as.data.frame(cbind(omega,keep_reps))
    names(tmp) <- c(paste0("omega",1:J),"iteration")
    return(tmp)
  },simplify = FALSE)
  mcmc <- do.call(rbind,mcmc)
  mcmc$chain <- factor(rep(1:chains,each=nrow(mcmc)/chains))
  mcmc <- mcmc %>% select(chain,iteration,everything())
  rm(counter)

  return(mcmc)
}
