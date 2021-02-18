#' Correlations between each pairwise feature and y
#'
#' @param pairwise_feature_mat a matrix with scores from feature pair comparison (rows) for each sample (cols)
#' @param y vector with quantitative outcome variable
#'
#' @return a sorted list of pairwise features based on their pearson correlation with y
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' y <- rnorm(C) # represents outcome variable
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' Is <- get_pairwise_rank_matrices(X)
#' pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
#' sorted_corrs <- get_sorted_corrs_pairwise_features(pairwise_feature_mat, y)
#'
get_sorted_corrs_pairwise_features <- function(pairwise_feature_mat, y){
  # Return a sorted list of pairwise features based on their pearson correlation with y 
  sorted_corrs <- sort(apply(pairwise_feature_mat, 2, function(x){
    stats::cor(x,y)
  }), decreasing = T)
  
  return(sorted_corrs)
}


#' Calculate weighted corr score
#'
#' Calculates a weighted correlation score that over weights higher correlations. When calling elders, this allows me to prioritize features that have are a part of a smaller number of highly correlative (with y) pairwise features as opposed to features that are a part of a larger number of relative less correlative pairwise features.
#'
#' @param corr a float between 0 and 1 representing the absolute value of a correlation coefficient (r)
#' @param max_score the maximum weight given to a correlation. Given that this weighted scoring can produce exponentially increasing weighted scores, setting this to 100 (representing a correlation of ~0.99) is reasonable.
#'
#' @return a float between 0 and 100, that represents a mapping of corr to a weighted score that overweights higher correlation values
#' @export
#'
#' @examples
#' corr <- 0.05
#' weighted_score <- elder_weighted_corr_score(corr)
#' 
elder_weighted_corr_score <- function(corr, max_score = 100){
  return(min(abs(1/log(corr)), max_score))
}

#' Calculate weighted feature correlation score
#'
#' For every feature, this function calculates a feature correlation score that is dervied from the sum of all weighted correlations (see weighting function) between y and pairwise features that said feature is a component of.
#'
#' @param sorted_corrs a vector of floats between -1 and 1 representing the the correlation coeficient between a pairwise feature and the outcome variable y. Output by get_sorted_corrs_pairwise_features().
#'
#' @return a numeric vector of length one. The value is the correlation of the elder with y when all pairwise features containing elder are accounted for (after weighting). The name of this numeric vector is the index of the elder in X as a character vector (i.e. string)
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' y <- rnorm(C) # represents outcome variable
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' Is <- get_pairwise_rank_matrices(X)
#' pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
#' sorted_corrs <- get_sorted_corrs_pairwise_features(pairwise_feature_mat, y)
#' elder <- get_elder(sorted_corrs)
#' 
get_elder <- function(sorted_corrs){
  
  # refernced below repeatedly
  sorted_corr_names <- names(sorted_corrs)
  
  # init emplty vector in which each original feature (i.e. 1/2 of pairwise feature) will be initialized to a value of zero
  feature_vec <- c() 
  for (feature in unique(unlist(sapply(sorted_corr_names, strsplit, "_")))){
    feature_vec[feature] <- 0
  }
  
  # for every feature, tally a sum of all weighted correlations between y and pairwise features that said feature is a component of. Note, this correlation will be weighted first before being added to the running sum
  for (i in 1:length(sorted_corrs)){
    feature1_feature2 <- sorted_corr_names[i]
    corr <- abs(as.numeric(sorted_corrs[i]))
    feature1 <- stringr::str_split_fixed(feature1_feature2, "_", n = 2)[1]
    feature2 <- stringr::str_split_fixed(feature1_feature2, "_", n = 2)[2]
    feature_vec[feature1] <- feature_vec[feature1] + elder_weighted_corr_score(corr)
    feature_vec[feature2] <- feature_vec[feature2] + elder_weighted_corr_score(corr)
  }
  
  # sort feature_vec
  sorted_feature_vec <- sort(feature_vec, decreasing = T)
  
  return(sorted_feature_vec[1]) # top feature is the elder
}

#' Determine elder and siblings that make cluster
#' 
#' This function identifies the best elder with which to define the best cluster for regression prediction. Such a cluster is made of pairwise features, the component features of which are termed elder and sibling, where the elder has the highest weighted sum of correlations with y across all pairwise features that it is a component of. Siblings are then the features that make up these pairwise features with elder.
#'
#' @param elder character vector (i.e. string) of the elder index I am up to
#' @param cluster_corr_prop kTSCR hyperparameter that determines what proportion oan elders overall correlation weight I want to capture in the included current cluster (which is made of elder sibling pairs). Default value = 1, which represents all elder sibling pairs.
#' @param sorted_corrs a vector of floats between -1 and 1 representing the the correlation coeficient between a pairwise feature and the outcome variable y. Output by get_sorted_corrs_pairwise_features().
#' @param elder_corr a length numeric vector (single value with a name) representing the sum weighted total correlation for the elder feature.
#'
#' @return a list of elder sibling pairs that represent the sibling pairs with elder that comprise the top correlations with y.
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
#' siblings <- get_siblings(sorted_corrs, elder, elder_corr)
#' 
get_siblings <- function(sorted_corrs, elder, elder_corr, cluster_corr_prop = 1){
  
  # now get siblings
  elder_corr_so_far <- 0 # init at zero - will record a sum of how much of the elder corr value has been accounted for as every sibling pair is found
  siblings <- c() # stores sibling as I find them
  for (i in 1:length(sorted_corrs)){
    
    # BREAK if elder_corr_so_far / elder_corr >  cluster_corr_prop
    # this is where hyperparameter cluster_corr_prop exerts its effect - it controls the number of siblings based on proportion fo total weighted elder correlation that has been captured so far
    if (elder_corr_so_far / elder_corr >= cluster_corr_prop){ 
      break
    }
    feature1_feature2 <- names(sorted_corrs)[i]
    feature1 <- stringr::str_split_fixed(feature1_feature2, "_", n = 2)[1]
    feature2 <- stringr::str_split_fixed(feature1_feature2, "_", n = 2)[2]
    # reverse the feature1_feature2 order if the correlation is negative
    corr <- as.numeric(sorted_corrs[i])
    if (corr < 0){
      feature1_feature2 <- paste(feature2, feature1, sep = "_")
    }
      
   if (feature1 == elder | feature2 == elder){
      siblings <- append(siblings, feature1_feature2)
      # add corr value to elder_val_build_up
      elder_corr_so_far <- elder_corr_so_far + elder_weighted_corr_score(abs(corr))
    }
  }
  
  return(siblings)
  
}

#' Get indices of sibling features
#'
#' Convert sibling feature names into indices from the original feature matrix X
#'
#' @param X input feature vector
#' @param siblings a vector of sibling names (each of which is a character vector, i.e. string)
#'
#' @return returns siblings as a numeric vector of indices rather than the character vector version of said indices
#' @export
#'
#' @examples
#' siblings <- c("V1_V3", "V20_V5", "V100_V110", "V101_V50")
#' sibling_indices <- get_sibling_indices(siblings)
#' 
#' 
get_sibling_indices <- function(X, siblings){
  
  if (is.null(rownames(X))){
    indices <- t(apply(do.call(rbind, strsplit(gsub("V", "", siblings), "_")), 1, as.numeric))
  } else{
    mapdf <- data.frame(feature_names = rownames(X), feature_indices = 1:nrow(X))
    indices <- matrix(mapdf$feature_indices[match(do.call(rbind, strsplit(tt, "_")),mapdf$feature_names)], ncol = 2, byrow = FALSE)
  }
  return(indices)
}