#' Combine per sample feature pair score matrices into a single non-redundant matrix (i.e. keep Xi - Xj, but not Xj - Xi, and have scores for every sample in a single matrix)  
#'
#' @param Is list of matrices, where each matrix represents a feature x feature score matrix for a single sample
#'
#' @return a matrix of feature pair score for each sample
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' Is <- get_pairwise_rank_matrices(X)
#' pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
#' 
make_feature_pair_score_matrix <- function(Is){
  num_features <- nrow(Is[[1]]) * (nrow(Is[[1]]) - 1) / 2 # size of non-redundant output matrix
  num_samples <- length(Is)
  pairwise_feature_mat <- matrix(rep(0, num_samples*num_features),
                  nrow = num_features) # init zero matrix of correct size
  
  new_var_names <- rep(NA, (nrow(Is[[1]])-1 * nrow(Is[[1]]) / 2)) # init new variable name vector of correct size
  feature_counter <- 1
  for (i in 1:(nrow(Is[[1]])-1)){
    for (j in (i+1):nrow(Is[[1]])){
      for (l in 1:length(Is)){
        pairwise_feature_mat[feature_counter,l] <- Is[[l]][j,i]
      }
      # set new var name
      new_var_name <- paste0(j, '_', i) # Delete me when refactoring
      # save new var name
      new_var_names[feature_counter] <- new_var_name
      # increment counter
      feature_counter <- feature_counter + 1
    }
  }
  
  # transpose matrix which sets variables as cols
  pairwise_feature_mat <- t(pairwise_feature_mat)
  # set column names of pairwise_feature_mat with new_var_names
  colnames(pairwise_feature_mat) <- new_var_names
  
  # if any column is entirely identical at every element,
  # add jitter to prevent issues when calculating correlations with y
  pairwise_feature_mat <- apply(pairwise_feature_mat, 2, function(col){
    if(all_equal(col)){
      return(col + stats::rnorm(length(col)))
    }else{
      return(col)
    }
  })
  
  return(pairwise_feature_mat)
}
