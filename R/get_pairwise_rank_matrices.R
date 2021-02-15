#' Calculate feature pair scores
#' 
#' This function compares the magnitude of the features in pairwise fashion (e.g Xi vs Xj). This is done across all (non-restricted) gene pairs per sample, and generates a list of num_samples matrices, where each matrix is an num_rows x num_rows matrix with either a) an indicator variable (1 if Xi > Xj) or b) a rank difference between features (rank(Xi) - rank(Xj))
#'
#' @param X a feature matrix with samples as cols and features as rows (NOTE: this is a transpose of how a feature matrix is frequently thought of)
#' @param restrict Defaults to FALSE. Otherwise, must be a vector of feature names to restrict analysis to. If provided, the output matrices will only contain pairwise feature scores for comparisons that always include at least one of the features from the restrict vector.In other words, the rows of the output matrices will be restricted to the features in restrict although the cols will still include all features This way, all restrict features are compared to all other features, which allows for the inclusion of non-restrict features in pairs, which may reflect important feature relationships. Otherwise, if restrict == false, all features will be considered.
#' @param use_diff_rank a logical. If FALSE, pairwise feature scores are indicators of whether Xi > Xj. If TRUE, pairwise feature scores are differences in rank
#'
#' @return a list of matrices, with each matrix containing a score reflecting a pairwise comparison between features, per sample (i.e. each matrix reflects pairwise feature scores per sample)
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' Is <- get_pairwise_rank_matrices(X)
#' 
get_pairwise_rank_matrices <- function(X, restrict = FALSE, use_diff_rank = TRUE){
  
  # first initializes zero matrix with dimensions num_rows x num_rows...used below to generate lower.triangle matrix
  S <- matrix(0, nrow = nrow(X), ncol = nrow(X)) # init r x r matrix
  
  # generate lower triangle matrix
  LT <- lower.tri(S) * 1 # lower triangle matrix - used below for rank-based feature comparison
  
  # rank features per sample - used to compare feature magnitudes in pairwise fashion
  ranks <- apply(X, 2, rank) # feature ranks per sample
  ranks <- split(ranks, rep(1:ncol(ranks), each = nrow(ranks))) # splits each col into a list
  
  if (isTRUE(use_diff_rank)){
    # use LT and ranks to generate a list of pairwise comparison matrices (one per sample)
    Is <- lapply(ranks, function(x){ # list of I matrices, where each I matrix is an R x R matrix with rank difference between features
      if (isFALSE(restrict)){
        sapply(x, function(b){x-b}) # x-b gives me the magnitude ordering I want
      } else{
        sapply(x, function(b){x-b})[rownames(X) %in% restrict, ]
      }
    })
  } else{
    # 0/1 versus rank difference - i.e. the old "indicator function" way
    Is <- lapply(ranks, function(x){ # list of I matrices, where each I matrix is an R x R matrix with 1 whereever feature has a bigger magnitude than other features
      if (isFALSE(restrict)){
        LT[x, x]
      } else{
        LT[x, x][rownames(X) %in% restrict, ]
      }
    })
  }
  
  # return
  return(Is)
}
