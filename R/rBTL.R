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
