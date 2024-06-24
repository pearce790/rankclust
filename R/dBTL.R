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
      terms <- unlist(lapply(1:Ri,function(r){theta[pi[r]]-logSumExp(theta[setdiff(pi[1:Ri],pi[0:(r-1)])])}))
      c(terms,rep(NA,J-Ri))
    }))
  }else{
    terms_log <- t(apply(Pi,1,function(pi){
      Ri <- length(na.exclude(pi))
      if(Ri == (J-1)){stop("Pi invalid: Cannot have Ri = J-1")}
      terms <- unlist(lapply(1:Ri,function(r){theta[pi[r]]-logSumExp(theta[setdiff(1:J,pi[0:(r-1)])])}))
      c(terms,rep(NA,J-Ri))
    }))
  }

  # return result
  if(log){return(sum(terms_log,na.rm=T))
  }else{return(exp(sum(terms_log,na.rm=T)))}
}
