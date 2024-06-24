mcmc_RCBTL <- function(Pi,J,a_gamma,b_gamma,lambda,nu0=NULL,num_iters=1000,nu_reps=2,groupwise=FALSE,
                       chains=4,burn_prop=0.5,thin=1,seed=NULL,normalize_omega=TRUE){

  if(!is.null(seed)){set.seed(seed)}
  counter <- 1
  mcmc <- replicate(n=chains,{
    print(paste0("Estimating chain ",counter," of ",chains,"."))
    counter <<- counter+1
    res <- fit_RCBTL(Pi=Pi,J=J,a_gamma=a_gamma,b_gamma=b_gamma,lambda=lambda,nu0=nu0,
                     num_iters=num_iters,nu_reps=nu_reps,groupwise=groupwise)
    nreps <- nrow(res$omega)
    keep_reps <- seq(ceiling(burn_prop*nreps)+1,nreps,by=thin)

    if(normalize_omega){
      omega <- t(apply(res$omega[keep_reps,,drop=FALSE],1,function(x){x/sum(x)}))
    }else{omega <- res$omega[keep_reps,,drop=FALSE]}
    tmp <- as.data.frame(cbind(omega,res$nu[keep_reps,,drop=FALSE],
                               res$g[keep_reps,,drop=FALSE],res$K[keep_reps],keep_reps))
    names(tmp) <- c(paste0("omega",1:J),paste0("nu",1:ncol(res$nu)),paste0("G",1:J),"K","iteration")
    return(tmp)
  },simplify = FALSE)
  mcmc <- do.call(rbind,mcmc)
  mcmc$chain <- factor(rep(1:chains,each=nrow(mcmc)/chains))
  mcmc <- mcmc %>% select(chain,iteration,K,everything())
  rm(counter)

  return(mcmc)
}
