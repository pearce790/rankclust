% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_partition.R
\name{sample_partition}
\alias{sample_partition}
\title{Partition sampling in Bayesian Rank-Clustered BTL models (internal use only)}
\usage{
sample_partition(
  Pi,
  I,
  J,
  nu,
  g,
  K,
  a_gamma,
  b_gamma,
  logprior_partition,
  b_g = 0.5,
  d_g = 0.5,
  groupwise = FALSE
)
}
\arguments{
\item{Pi}{A matrix of preference orderings ("rankings"), such that the (i,j) entry is the index of the jth-most preferred object according to judge i. If \code{groupwise=TRUE}, then the index corresponds to the jth-most preferred object among those in row i; if \code{groupwise=FALSE}, it is assumed that all unranked objects (if any) are less preferred than those which are ranked.}

\item{I}{A numeric indicating the number of rows in Pi}

\item{J}{A numeric indicating the total number of objects being compared.}

\item{nu}{A vector indicating current values for nu in the Gibbs sampler.}

\item{g}{A vector indicating current values for g in the Gibbs sampler.}

\item{K}{A vector indicating current values for K in the Gibbs sampler.}

\item{a_gamma}{A numeric for the first hyperparameter (shape) in a Gamma prior on each worth parameter.}

\item{b_gamma}{A numeric for the second hyperparameter (rate) in a Gamma prior on each worth parameter.}

\item{logprior_partition}{A J-vector indicating the log unnormalized probability for each possible partition, based on its number of non-empty clusters.}

\item{b_g}{The probability of "birth"ing a new partition cluster, if possible.}

\item{d_g}{The probability of "death"ing an existing partition cluster, if possible.}

\item{groupwise}{A boolean to indicate whether the observed rankings are complete/partial rankings (\code{FALSE}; default) or groupwise comparisons (\code{TRUE}).}
}
\value{
A list containing updated values for g, nu, and K.
}
\description{
This function implements a reversible jump MCMC procedure for updating the parameter partition in Bayesian Rank-Clustered BTL models. For internal use only.
}
