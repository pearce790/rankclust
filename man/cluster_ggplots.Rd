% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_ggplots.R
\name{cluster_ggplots}
\alias{cluster_ggplots}
\title{Rank-Cluster ggplots for Rank-Clustered BTL models.}
\usage{
cluster_ggplots(mcmc, object_names = NULL)
}
\arguments{
\item{mcmc}{The matrix output of the \code{mcmc_RCBTL} function.}

\item{object_names}{An optional vector of object names, in the same index order as the observed data.}
}
\value{
A ggplot(s) of the pairwise rank-clustering posterior probabilities.
}
\description{
This function creates a rank-clustering probability matrix ggplot using MCMC chains from a Rank-Clustered BTL distribution.
}
\examples{
Pi <- rBTL(I=30,omega=c(4^2,4^2,4^1,4^1,4^0))
mcmc <- mcmc_RCBTL(Pi=Pi,J=5,a_gamma=5,b_gamma=3,lambda=2,num_iters=200,nu_reps=1,chains=2,seed=1)
cluster_ggplots(mcmc,object_names=paste0("Object ",1:5))
}
