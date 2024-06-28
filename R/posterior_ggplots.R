#' Posterior ggplots for worth parameters in (Rank-Clustered) BTL models.
#'
#' This function creates side-by-side violin plots for omega (object-specific worth parameters) using MCMC chains from a BTL or Rank-Clustered BTL distribution.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @param mcmc The matrix output of the \code{mcmc_BTL} or \code{mcmc_RCBTL} functions.
#' @param object_names An optional vector of object names, in the same index order as the observed data.
#'
#' @return A ggplot(s) of side-by-side violin plots for each omega_j.
#'
#' @examples
#' Pi <- rBTL(I=30,omega=c(4^2,4^2,4^1,4^1,4^0))
#' mcmc <- mcmc_BTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,lambda=2,num_iters=200,nu_reps=1,chains=2,seed=1)
#' posterior_ggplots(mcmc,object_names=paste0("Object ",1:5))
#' @export
posterior_ggplots <- function(mcmc,object_names=NULL){
  names_mcmc <- setdiff(names(mcmc),c("chain","iteration"))

  omega_vars <- grep("omega",names_mcmc)
  J <- length(omega_vars)
  if(is.null(object_names)){
    posteriorOmega <- ggplot(reshape2::melt(mcmc[,c("chain","iteration",names_mcmc[omega_vars])],id.vars=c(1,2)),
                         aes(x=factor(variable,levels=paste0("omega",1:J),labels=1:J),y=value))+
      geom_violin()+theme_bw()+
      labs(x=element_blank(),y=expression(omega))
  }else{
    posteriorOmega <- ggplot(reshape2::melt(mcmc[,c("chain","iteration",names_mcmc[omega_vars])],id.vars=c(1,2)),
                             aes(x=factor(variable,levels=paste0("omega",1:J),labels=object_names),y=value))+
      geom_violin()+theme_bw()+
      labs(x=element_blank(),y=expression(omega))
  }
  return(posteriorOmega)
}
