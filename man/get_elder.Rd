% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_siblings.R
\name{get_elder}
\alias{get_elder}
\title{Calculate weighted feature correlation score}
\usage{
get_elder(sorted_corrs)
}
\arguments{
\item{sorted_corrs}{a vector of floats between -1 and 1 representing the the correlation coeficient between a pairwise feature and the outcome variable y. Output by get_sorted_corrs_pairwise_features().}
}
\value{
a numeric vector of length one. The value is the correlation of the elder with y when all pairwise features containing elder are accounted for (after weighting). The name of this numeric vector is the index of the elder in X as a character vector (i.e. string)
}
\description{
For every feature, this function calculates a feature correlation score that is dervied from the sum of all weighted correlations (see weighting function) between y and pairwise features that said feature is a component of.
}
\examples{
C <- 100  # represents samples
R <- 200 # represents features
y <- rnorm(C) # represents outcome variable
X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
Is <- get_pairwise_rank_matrices(X)
pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
sorted_corrs <- get_sorted_corrs_pairwise_features(pairwise_feature_mat, y)
elder <- get_elder(sorted_corrs)

}
