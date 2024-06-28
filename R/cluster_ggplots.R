#' Rank-Cluster ggplots for Rank-Clustered BTL models.
#'
#' This function creates a rank-clustering probability matrix ggplot using MCMC chains from a Rank-Clustered BTL distribution.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @param mcmc The matrix output of the \code{mcmc_RCBTL} function.
#' @param object_names An optional vector of object names, in the same index order as the observed data.
#'
#' @return A ggplot(s) of the pairwise rank-clustering posterior probabilities.
#'
#' @examples
#' Pi <- rBTL(I=30,omega=c(4^2,4^2,4^1,4^1,4^0))
#' mcmc <- mcmc_RCBTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,lambda=2,num_iters=200,nu_reps=1,chains=2,seed=1)
#' cluster_ggplots(mcmc,object_names=paste0("Object ",1:5))
#' @export
cluster_ggplots <- function(mcmc,object_names=NULL){
  names_mcmc <- setdiff(names(mcmc),c("chain","iteration"))
  if(length(grep("G",names_mcmc))==0){
    stop("Incorrect input. Must be the return object of the `mcmc_RCBTL` function.")
  }

  posterior_median_order <- rev(order(apply(mcmc[,names_mcmc[grep("omega",names_mcmc)]],2,median)))
  J <- length(posterior_median_order)
  pairs <- expand.grid(i=1:J,j=1:J,prob=0)
  pairs <- pairs[pairs$i < pairs$j,]
  pairs$prob <- apply(pairs,1,function(pair){
    i <- pair[1]; j <- pair[2]
    mean(apply(mcmc[,paste0("G",c(i,j))],1,function(x){x[1]==x[2]}))
  })
  pairs2 <- pairs[,c(2,1,3)]; names(pairs2) <- c("i","j","prob")
  pairs3 <- data.frame(i=1:J,j=1:J,prob=1)
  pairs <- rbind(pairs,pairs2,pairs3)

  if(is.null(object_names)){
    clusters <- ggplot(pairs,aes(x=factor(i,levels=posterior_median_order),
                                 y=factor(j,levels=rev(posterior_median_order)),
                                 fill=prob))+
      geom_tile()+scale_fill_gradient(low="white",high="black",breaks=c(0,.5,1))+
      theme_bw()+labs(x=element_blank(),y=element_blank(),fill="Probability")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }else{
    clusters <- ggplot(pairs,aes(x=factor(i,levels=posterior_median_order,labels=object_names[posterior_median_order]),
                                 y=factor(j,levels=rev(posterior_median_order),labels=object_names[rev(posterior_median_order)]),
                                 fill=prob))+
      geom_tile()+scale_fill_gradient(low="white",high="black",breaks=c(0,.5,1))+
      theme_bw()+labs(x=element_blank(),y=element_blank(),fill="Probability")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(clusters)
}
