#' Implements kTSCR algorithm
#' 
#' This algorithm find the top clusters for predicting a quantitative outcome variable (y) using input features (X). Clusters are defined by a singular feature, called an elder, and all other features that optimize prediction when compared with the elder (either as an indicator function, or as difference of rank). As such, the elder-sibling pairs represent transformed features (termed pairwise features), and all elder-sibling pairs comprise the cluster.
#'
#' @param y vector with quantitative outcome variable
#' @param X a feature matrix with samples as cols and features as rows (NOTE: this is a transpose of how a feature matrix is frequently thought of)
#' @param Verbose a logical. Whether tobe verbose (TRUE) or not. Default is TRUE.
#' @param restrict Defaults to FALSE. Otherwise, must be a vector of feature names to restrict analysis to. If provided, the output matrices will only contain pairwise feature scores for comparisons that always include at least one of the features from the restrict vector.In other words, the rows of the output matrices will be restricted to the features in restrict although the cols will still include all features This way, all restrict features are compared to all other features, which allows for the inclusion of non-restrict features in pairs, which may reflect important feature relationships. Otherwise, if restrict == false, all features will be considered.
#' @param rank a logical. Whether the rank of the outcome variable should be used for prediction. Default is FALSE. Recommend against using TRUE for this parameter (currently used for testing).
#' @param standardize_features a logical. Whether features in X should be standardized. Default is TRUE.
#' @param cluster_corr_prop what proportion of the maximum (weighted) cluster correlation with y should be reflected by the chosen siblings. A hyperparameter. Default is 1 (meaning include all elder-sibling pairs in cluster)
#' @param ct correlation threshold determined how much a new cluster must improve the current correlation with y in order to be added as a top cluster. A hyperparameter. Default is 1 (meaning any improvement is sufficient to add the next cluster within the greedy framework)
#' @param use_diff_rank a logical. If true feature pairs are scored based on the difference in their per sample rank. Otherwise, pairwise scores are the indicator (I) of whether Xi>Xj
#'
#' @return a list with entries y (outcome variable), K (calculated score used for prediction), Correlation (cor(K, y)), ElderIndices (indices of elders in X), SiblingIndices (indices of siblings in X), Elders (elder variable names), siblings (sibling variable names)
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' y <- rnorm(C) # represents outcome variable
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' res <- get_top_clusters(y, X)
#' 
get_top_clusters <- function(y,
                             X,
                             Verbose = FALSE,
                             restrict = FALSE,
                             rank = FALSE,
                             standardize_features = TRUE,
                             cluster_corr_prop = 1,
                             ct = 1.0, # 1.01 would mean that correlation must improve by more than 1% to keep searching for elders/clusters during greedy optimization
                             use_diff_rank = T
){
  
  # standardize feature matrix if standardize_features
  # argument == T
  # i.e. scale and center each feature
  if (isTRUE(standardize_features)){
    if (Verbose){
      print('Standardizing features', quote = F)
    }
    X <- t(apply(X, 1, scale, center = T, scale = T)) # must transpose back to original dimensions due to how apply returns matrix object
  }
  
  # use rank of outcome if rank argument specified
  if (isTRUE(rank)){
    y <- rank(y)
  }
  
  # require that length(y) == ncol(X)
  stopifnot(length(y) == ncol(X))
  
  # INITIALIZE VARIABLES
 
  if (Verbose){
    
    print('Initializing variables', quote = F)
  }
  
  K <- rep(0, ncol(X))  # init K <- 0 for every sample (stores running K score)
  G <- list() # init list of top hit features (termed elders)
  B <- list() # init list to hold list of siblings for each elder
  k <- 1 # little k - pair you're up to - max K must be <= k each iteration
  last_corr <- 0
  
  # PERFORM INITIAL CALCULATIONS
  
  # generate list of  pairwise comparison matrices
  if (Verbose){
    print('Calculating pairwise indicator matrices', quote = F)
  }
  
  Is <- get_pairwise_rank_matrices(X, restrict, use_diff_rank)
  pairwise_feature_mat <- make_feature_pair_score_matrix(Is)
  
  
  # WHILE LOOP
  
  # use a while loop to identify pairs that optimize correlation
  
  if (Verbose){
    print('Starting kTSCR iteration', quote = F)
  }
  
  while(TRUE){
    
    if (Verbose){
      print(paste0("    Iteration ", k), quote = F)
    }
    
    # if all features have already been found to contain a brother, break
    if (dim(pairwise_feature_mat)[2] == 0){
      if (Verbose){
        print(paste0("Converged at k=", k,": ALL siblings"), quote = F)
      }
      break
    }
    
    if (dim(pairwise_feature_mat)[2] == 1){ # TODO: fix this - one way is to just choose the pair and see if adding it improves the correlation
    }
    
    # FIND ELDER AND SIBLINGS
    
    # first calculate the sorted corrs
    sorted_corrs <- get_sorted_corrs_pairwise_features(pairwise_feature_mat, y)
    # identify ELDER using correlation-based feature scores - this is done each iteration using the updated pairwise_feature_mat
    elder_corr <- get_elder(sorted_corrs)
    elder <- names(elder_corr)
    # Use elder to identify SIBLINGS (elder + siblings define the cluster)
    siblings <- get_siblings(sorted_corrs, elder, elder_corr, cluster_corr_prop = cluster_corr_prop)
    
    # Check to see if any siblings were found.
    # If they weren't, then BREAK the while loop and DO NOT
    # save the elder or any updates to any variables
    # i.e. the greedy optimization has finished
    if (is.null(siblings)){ # zero length vector (c()) evaluates to NULL
      if (Verbose){
        print(paste0("Converged at k=", k,": NO siblings"), quote = F)
      }
      break
    }
    
    # if you don't break...
    # convert sibling pairs into indices for use in calculating K score with them
    
    siblings <- split_sibling_indices(siblings)
    
    # Test for another WHILE LOOP BREAK CONDITION - CORRELATION CONVERGENCE
    new_K <- K + calc_K_from_pairs(siblings, X, restrict = restrict)
    current_corr <- stats::cor(new_K, y)
    if (current_corr <= (last_corr*ct)){ # ct - convergenece tolerance parameter (see above)
      if (Verbose){
        print(paste0("Converged at k=", k,": CORRELATION DOES NOT SIGNIFICANTLY IMPROVE"))
      }
      break
    }
    
    # print correlation thus far
    if (Verbose){
      print(paste0("      Elder: ", elder, " Correlation: ", current_corr), quote = F)
    }
    
    # UPDATE PARAMETER PRIOR TO NEXT ITERATION
    
    K <- new_K # update K score
    G[[k]] <- elder # save chosen Elder
    B[[k]] <- siblings
    
    # update k for next iteration (remember, in this implementation, k is one ahead whenever I break the while loop)
    k <- k + 1
    
    # break if k is already equal to the number of total features
    # (this should very rarely happen in practice), as there would be
    # no more elders to find in the next iteration
    if (k > nrow(X)){
      if(Verbose){
        print(paste0("Converged at k=", k,": ALL FEATURES CHOSEN AS ELDERS"))
      }
      break
    }
    
    # update last_corr
    last_corr <- current_corr
    
    # update pairwise_feature_mat
    pairwise_feature_mat <- update_pairwise_feature_mat(pairwise_feature_mat, elder)
  }
  
  # OUTPUT RELEVANT INFORMATION
  
  # Now that you are outside the while loop,
  # return relevant information as a list
  # Note: look further into customized object to return data in
  
  # store indices of elders and brother
  elder_indices <- as.numeric(unlist(G))
  sibling_indices <- do.call(rbind, B)
  
  if (!is.null(rownames(X))){
    elders <- sapply(G, function(v){ # return vector of elders
      i <- as.numeric(v)
      rownames(X)[i]
    })
    siblings <- do.call(rbind,
                        lapply(B, function(v){
                          i <- as.numeric(v)
                          matrix(rownames(X)[i], byrow = F, nrow = nrow(v))
                        }))
  } else{
    elders <- NULL
    siblings <- NULL
  }
  
  # list to output
  res_list <- list(y = y,
                  K = K,
                  ApparentCorrelation = stats::cor(K, y),
                  ElderIndices = elder_indices,
                  SiblingIndices = sibling_indices,
                  Elders = elders,
                  Siblings = siblings
  )
  
  # return
  return(res_list)
  
  if (Verbose){
    print('DONE', quote = F)
  }
}
