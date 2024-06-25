#' Random data generation from standard BTL models
#'
#' This function randomly draws data from a Bradley-Terry-Luce (BTL) model with object-specific worth parameters, \code{omega}. Currently, only complete or partial rankings may be drawn, although those may be subsetted to simulate groupwise comparison data as a consequence of Luce's Axiom of Choice.
#'
#' @param I A numeric indicating the number of observations to draw.
#' @param omega A vector of non-negative object worth parameters.
#' @param R A numeric indicating the length of each ranking. Default to \code{NULL}, indicating to draw complete rankings.
#'
#' @return A \code{I}x\code{R} matrix of rankings.
#'
#' @examples
#' set.seed(1)
#' rBTL(I=5,omega=c(10,5,1),R=3)
#' @export
rBTL <- function(I,omega,R=NULL){

  if((I %% 1)!=0 | I<1){stop("I must be a positive integer.")}
  if(any(omega<0)){stop("All entries in omega must be non-negative.")}
  J <- length(omega)
  Pi <- t(replicate(I,{sample(1:J,J,replace=FALSE,prob=omega)}))

  if(!is.null(R)){
    if(any(R%%1 != 0) | any(R<1) | any(R>J)){stop("R must be an integer between 1 and J")}
    if(any(R==(J-1))){stop("R cannot equal J-1; this is equivalent to R=J")}
    if(length(R)==1){Pi <- Pi[,1:R,drop=FALSE]
    }else if(length(R)==I){
      for(i in 1:I){if(R[i]<J){Pi[i,(R[i]+1):J] <- NA}}
      Pi <- Pi[,1:max(R),drop=FALSE]
    }else{stop("R must be of length 1 or I")}
  }

  if(I==1){return(matrix(Pi,nrow=1))
  }else{return(Pi)}
}
