#' Calculate K score from pairs
#'
#' calculate K score, which should correlate with y and therefore serve as the basis for the regression based prediction. Use the input feature matrix (X) and a set of identified of identified (elder-sibling) pairs to calculate K
#'
#' @param pairs Pairwise feature pairs (Elder-sibling) identified during optimization procedure (get_siblings())
#' @param X Input feature matrix
#' @param restrict Defaults to FALSE. Otherwise, must be a vector of feature names to restrict analysis to. If provided, the output matrices will only contain pairwise feature scores for comparisons that always include at least one of the features from the restrict vector.In other words, the rows of the output matrices will be restricted to the features in restrict although the cols will still include all features This way, all restrict features are compared to all other features, which allows for the inclusion of non-restrict features in pairs, which may reflect important feature relationships. Otherwise, if restrict == false, all features will be considered. Used below by the get_pairwise_rank_matrices() function, which while not technically necessary, is effectively necessary frmo a performance standpoint.
#'
#' @return a numeric vector with a score (K) for each sample (the correlation based score used for prediction that kTSCR calculates)
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' y <- rnorm(C) # represents outcome variable
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' pairs <- rbind(c(1,2),c(4,6), c(10,5))
#' K <- calc_K_from_pairs(pairs, X)
#' 
calc_K_from_pairs <- function(pairs, X, restrict = FALSE){
  
  # calculate K score given input feature matrix (X) and a set of
  # identified pairs
  # X is feature matrix to use for calculating K Score
  
  # convert X to Is
  Is <- get_pairwise_rank_matrices(X, restrict)
  
  K <- unlist(lapply(Is, function(j){ # iterate over matrices in Is list
    # j is a matrix from Is (represents a samples rank diff matrix)
    diffs <- apply(pairs, 1, function(i){ # iterates over feature pair indices
      return(j[i[1], i[2]]) # the way the pairs should be ordered
    })
    return(sum(diffs))
  }))
  
  return(K)
  
}
