% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_siblings.R
\name{get_siblings}
\alias{get_siblings}
\title{Determine elder and siblings that make cluster}
\usage{
get_siblings(sorted_corrs, elder, elder_corr, cluster_corr_prop = 1)
}
\arguments{
\item{sorted_corrs}{a vector of floats between -1 and 1 representing the the correlation coeficient between a pairwise feature and the outcome variable y. Output by get_sorted_corrs_pairwise_features().}

\item{elder}{character vector (i.e. string) of the elder index I am up to}

\item{elder_corr}{a length numeric vector (single value with a name) representing the sum weighted total correlation for the elder feature.}

\item{cluster_corr_prop}{kTSCR hyperparameter that determines what proportion oan elders overall correlation weight I want to capture in the included current cluster (which is made of elder sibling pairs). Default value = 1, which represents all elder sibling pairs.}
}
\value{
a list of elder sibling pairs that represent the sibling pairs with elder that comprise the top correlations with y.
}
\description{
This function identifies the best elder with which to define the best cluster for regression prediction. Such a cluster is made of pairwise features, the component features of which are termed elder and sibling, where the elder has the highest weighted sum of correlations with y across all pairwise features that it is a component of. Siblings are then the features that make up these pairwise features with elder.
}
\examples{
C <- 100  # represents samples
R <- 200 # represents features
y <- stats::rnorm(C) # represents outcome variable
X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
Is <- get_pairwise_rank_matrices(X)
pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
sorted_corrs <- get_sorted_corrs_pairwise_features(pairwise_feature_mat, y)
elder_corr <- get_elder(sorted_corrs)
elder <- names(elder_corr)
siblings <- get_siblings(sorted_corrs, elder, elder_corr)

}
