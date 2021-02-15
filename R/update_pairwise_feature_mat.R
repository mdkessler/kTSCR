#' Remove elder from pairwise feature matrix
#'
#' Remove all columns from pairwise_feature_mat that have elder as one of the features in each pairwise feature
#'
#' @param pairwise_feature_mat a matrix of feature pair score for each sample
#' @param elder character vector (i.e. string) representing the name of the elder feature
#'
#' @return returns pairwise feature matrix but without the columns that had pairwise features where one feature was the elder feature
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' y <- stats::rnorm(C) # represents outcome variable
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' Is <- get_pairwise_rank_matrices(X)
#' pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
#' sorted_corrs <- get_sorted_corrs_pairwise_features(pairwise_feature_mat, y)
#' elder_corr <- get_elder(sorted_corrs)
#' elder <- names(elder_corr)
#' pairwise_feature_mat <- update_pairwise_feature_mat(pairwise_feature_mat, elder)
#' 
update_pairwise_feature_mat <- function(pairwise_feature_mat, elder){
  # Remove all columns from pairwise_feature_mat that have elder as one of the
  # indicator comparison pairs
  
  exclude <- which(sapply(1:ncol(pairwise_feature_mat), function(i){
    elder %in% unlist(strsplit(colnames(pairwise_feature_mat)[i], "_"))
  }))
  
  # deal with the special case where there is only one indicator feature column left in the reduced (i.e. updated) matrix
  if (ncol(pairwise_feature_mat) - length(exclude) == 1){
    ret <- as.matrix(pairwise_feature_mat[, -exclude]) # as matrix is necessary to explictly convert 1 column matrix from vector back to matrix
    colnames(ret) <- colnames(pairwise_feature_mat)[!(1:ncol(pairwise_feature_mat)) %in% exclude]
    return(ret)
  }
  # else
  return(pairwise_feature_mat[, -exclude]) 
}
