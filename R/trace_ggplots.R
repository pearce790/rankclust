#' Trace ggplots for (Rank-Clustered) BTL models.
#'
#' This function creates a trace plot for K (the number of rank-clusters, if applicable) and omega (object-level worth parameters) using MCMC chains from a BTL or Rank-Clustered BTL distribution.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @param mcmc The matrix output of the \code{mcmc_BTL} or \code{mcmc_RCBTL} functions.
#' @param object_names An optional vector of object names, in the same index order as the observed data.
#'
#' @return A list containing ggplot(s) containing trace plots of K (if applicable) and omega.
#'
#' @examples
#' Pi <- rBTL(I=30,omega=c(4^2,4^2,4^1,4^1,4^0))
#' mcmc <- mcmc_RCBTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,lambda=2,num_iters=200,nu_reps=1,chains=2,seed=1)
#' trace_ggplots(mcmc,object_names=paste0("Object ",1:5))
#' @export
trace_ggplots <- function(mcmc,object_names=NULL){
  names_mcmc <- setdiff(names(mcmc),c("chain","iteration"))

  K_vars <- grep("K",names_mcmc)
  if(length(K_vars)==0){traceK <- NULL
  }else{
    traceK <- ggplot(mcmc[,c("chain","iteration",names_mcmc[K_vars])],
                     aes(x=iteration,y=K,color=chain))+
      geom_line()+theme_bw()+
      labs(x="Iteration",y="K",color="Chain")
  }

  omega_vars <- grep("omega",names_mcmc)
  J <- length(omega_vars)
  if(is.null(object_names)){
    traceOmega <- ggplot(reshape2::melt(mcmc[,c("chain","iteration",names_mcmc[omega_vars])],id.vars=c(1,2)),
                         aes(x=iteration,y=value,color=chain))+
      facet_wrap(~factor(variable,levels=paste0("omega",1:J),labels=1:J))+
      geom_line()+theme_bw()+
      labs(x="Iteration",y=expression(omega),color="Chain")
  }else{
    traceOmega <- ggplot(reshape2::melt(mcmc[,c("chain","iteration",names_mcmc[omega_vars])],id.vars=c(1,2)),
                         aes(x=iteration,y=value,color=chain))+
      facet_wrap(~factor(variable,levels=paste0("omega",1:J),labels=object_names))+
      geom_line()+theme_bw()+
      labs(x="Iteration",y=expression(omega),color="Chain")
  }

  result <- list(traceK=traceK,traceOmega=traceOmega)
  return(result)
}
