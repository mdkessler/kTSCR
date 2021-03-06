% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_cross_validation.R
\name{run_cross_validation}
\alias{run_cross_validation}
\title{Run k fold cross validation}
\usage{
run_cross_validation(
  y,
  X,
  Verbose.pass = FALSE,
  restrict = FALSE,
  rank = FALSE,
  Verbose = TRUE,
  standardize_features = TRUE,
  cluster_corr_prop = 1,
  ct = 1,
  sibling_prune = 0.1,
  k = 5,
  condensed_output = TRUE
)
}
\arguments{
\item{y}{outcome variable}

\item{X}{inut feature matrix}

\item{Verbose.pass}{logical as to whether the kTSCR procedure should be verbose (i.e. should run_cross_validation pass 'verbose=TRUE' to get_top_clusters()) )}

\item{restrict}{a list of colnames of X by which to restrict the analysis}

\item{rank}{logical as to whether to use rank of outcome}

\item{Verbose}{logical as to whether to be verbose}

\item{standardize_features}{logical as to whether to standardize all features of X}

\item{cluster_corr_prop}{what proportion of the maximum (weighted) cluster correlation with y should be reflected by the chosen siblings. A hyperparameter. Default is 1 (meaning include all elder-sibling pairs in cluster)}

\item{ct}{correlation threshold determined how much a new cluster must improve the current correlation with y in order to be added as a top cluster. A hyperparameter. Default is 1 (meaning any improvement is sufficient to add the next cluster within the greedy framework)}

\item{sibling_prune}{numeric between 0-1 that sets the threshold for how close apparent correlation and test correlation must be for a k-cv iteration to contribute its siblings to the final chosen siblings. In other words, a lower number is more stringent, since it means the overfitting had to be really low in a k-cv iteration for it to contribute to the final sibling output.}

\item{k}{the k parameter in k fold cross validation (i.e. train/test partitions). Default is 5}

\item{condensed_output}{return output that is condensed and summarized across k_cv iterations, specifically with regard to feature importance}
}
\value{
returns the list given by get_top_clusters for each n fold k cv run and includes the test correlation and train/test splits from each iteration
}
\description{
This function runs k fold cross validation by splitting input data into k partitions and holding out each partition as the test set in k different learning iterations.
}
\examples{
C <- 100  # represents samples
R <- 200 # represents features
y <- rnorm(C) # represents outcome variable
X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
cv_res <- run_cross_validation

}
