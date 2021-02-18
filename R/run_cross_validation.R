#' Calculate apparent correlation minus test correlation
#'
#' Calculate app_cor - test_cor as a measure of overfitting at each k-cv iteration
#'
#' @param k_cv_res result from run_cross_validation (called with condensed_output = TRUE)
#' @param k k number of cross validation iterations, passed from run_cross_validation function
#'
#' @return a numeric vector with app_cor - test_cor
#' @export
#'
#' @examples
#' 
get_app_minus_test <- function(k_cv_res, k){
  app_corr_vec <- unlist(k_cv_res['app_cor', 1:k])
  test_corr_vec <- unlist(k_cv_res['test_cor', 1:k])
  app_minus_test <- app_corr_vec - test_corr_vec
  
  return(app_minus_test)
}

#' Weight feature importance sum on the basis of overfitting
#'
#' When calculating feature importance by tallying the number of times a feature is in a sibling within and then across k-cv iterations, weight the tally by the overfitting seen in the iteration, as measured by app_cor - test_cor
#'
#' @param k_cv_res result from run_cross_validation (called with condensed_output = TRUE)
#' @param feature_importance_mat feature importance values from each kTSCR iteration derived from how many times a feature was seen in a sibling
#' @param k k number of cross validation iterations, passed from run_cross_validation function
#' @param penalty a penalty used in the weighting. A (potential) hyperparameter. Default is 5.
#' 
#' @return numeric vector of weighted feature importance sums
#' @export
#'
#' @examples
#' 
app_minus_test_thresh_weighted_sum <- function(k_cv_res, feature_importance_mat, k, penalty = 5){
  app_minus_test <- get_app_minus_test(k_cv_res, k)
  app_minus_test_weights <- 1 / (1 + app_minus_test)**penalty # this ensure that overfit iterations are down-weighted and underfit iterations are upweighted
  # now return a linear combination of columns in feature_importance_mat weighted by app_minus_test_weights
  return(feature_importance_mat %*% app_minus_test_weights)
}  

#' Condense k cross validation output
#'
#' Condense and summarize output from k cross validation iterations. Specifically, this returns the union of elders across k-cv iterations, and the total per feature of how many sibling pairs it was in across all k-cv iterations
#'
#' @param X input feature matrix
#' @param k_cv_res k cross validation result object from function run_cross_validation()
#' @param k k number of cross validation iterations, passed from run_cross_validation function
#' @param app_minus_test_thresh a threshold level for app_cor - test_cor...k-cv iterations that have app_cor - test_cor <= app_minus_test_thresh (i.e. low overfitting) have their siblings added to the consensus siblings ultimately used for prediction. A hyperparameter. Defaults to 0.10
#' @param weight_sum_by_app_minus_thresh A logical. Indicates whether the sums across k-cv iterations used for feature importance should be weighted by app_cor - test_cor (i.e. higher weight to k-cv iterations that overfit less). Default is true.
#'
#' @return list of four vectors: app_cor, test_cor, elders, feature_importance
#' @export
#'
#' @examples
#' 
condense_k_cv_output <- function(X, k_cv_res, k, app_minus_test_thresh = 0.10, weight_sum_by_app_minus_thresh = TRUE){
  
  # first condense elders by getting the union of elders across k_cv iterations
  elders_vec <- unique(unlist(k_cv_res['elders', 1:k]))
  
  # now condense feature importance scores but summing them across k_cv iterations
  feature_importance_list <- k_cv_res['feature_importance', 1:k]
  feature_names <- names(feature_importance_list[[1]])
  feature_importance_mat <- matrix(unlist(feature_importance_list), ncol=length(feature_importance_list), byrow = FALSE)
  if (isTRUE(weight_sum_by_app_minus_thresh)){
    feature_importance_vec <- app_minus_test_thresh_weighted_sum(k_cv_res, feature_importance_mat, k)
  }else{
    feature_importance_vec <- apply(feature_importance_mat, 1, sum)
  }
  rownames(feature_importance_vec) <- feature_names
  
  # use feature importance scores to choose top features - this will be used below to filter siblings
  top_features <- choose_top_features(feature_importance_vec)
    
  # now condense sibling pairs
  
  # first, get the intersect of siblings across k-cv iterations
  siblings_list <- lapply(k_cv_res['siblings', 1:k], function(x){ # format feature names and paste features together per pairwise feature
    apply(x, 1, function(x){paste0(x, collapse = "_")})
  })
  siblings_vec <- Reduce(intersect, siblings_list) # get intersect
  
  # now add sibling from any k-cv iteration where app_cor - test_cor <= app_minus_test_thresh
  app_minus_test <- get_app_minus_test(k_cv_res, k)
  app_minus_test_pass_thresh <- app_minus_test <= app_minus_test_thresh
  if (sum(app_minus_test_pass_thresh) > 0){
    for (i in which(app_minus_test_pass_thresh)){
      # first, only keep sibling if one feature is an elder
      siblings_with_top_features <- which(
        k_cv_res['siblings', i][[1]][,1] %in% top_features
          | # or
        k_cv_res['siblings', i][[1]][,2] %in% top_features
      )
      
      siblings_filtered <- k_cv_res['siblings', i][[1]][siblings_with_top_features, ]
      # now format the names
      siblings_to_add <- apply(siblings_filtered, 1, function(x){paste0(x, collapse = "_")})
      siblings_vec <- union(siblings_vec, siblings_to_add)
    }
  } else{
    for (i in which(k_cv_res$test_corr >= 0.2)){
      # first, only keep sibling if one feature is an elder
      siblings_with_top_features <- which(
        k_cv_res['siblings', i][[1]][,1] %in% top_features
        | # or
          k_cv_res['siblings', i][[1]][,2] %in% top_features
      )
      
      siblings_filtered <- k_cv_res['siblings', i][[1]][siblings_with_top_features, ]
      # now format the names
      siblings_to_add <- apply(siblings_filtered, 1, function(x){paste0(x, collapse = "_")})
      siblings_vec <- union(siblings_vec, siblings_to_add)
    }
  }
  # convert siblings_vec back into indices matrix
  siblings_mat <- split_sibling_names(siblings_vec)
  siblings_indices_mat <- convert_sibling_names_to_indices(X, siblings_mat)
  
  # output vecs as a list
  app_corr_vec <- unlist(k_cv_res['app_cor', 1:k])
  test_corr_vec <- unlist(k_cv_res['test_cor', 1:k])
  return(list(app_corr = app_corr_vec, test_corr = test_corr_vec, elders = elders_vec, siblings = siblings_mat, siblings_indices = siblings_indices_mat, feature_importance = feature_importance_vec))
}


#' Run k fold cross validation
#'
#' This function runs k fold cross validation by splitting input data into k partitions and holding out each partition as the test set in k different learning iterations. 
#'
#' @param y outcome variable
#' @param X inut feature matrix
#' @param Verbose.pass logical as to whether the kTSCR procedure should be verbose (i.e. should run_cross_validation pass 'verbose=TRUE' to get_top_clusters()) )
#' @param restrict a list of colnames of X by which to restrict the analysis
#' @param rank logical as to whether to use rank of outcome
#' @param Verbose logical as to whether to be verbose
#' @param standardize_features logical as to whether to standardize all features of X
#' @param cluster_corr_prop what proportion of the maximum (weighted) cluster correlation with y should be reflected by the chosen siblings. A hyperparameter. Default is 1 (meaning include all elder-sibling pairs in cluster)
#' @param ct correlation threshold determined how much a new cluster must improve the current correlation with y in order to be added as a top cluster. A hyperparameter. Default is 1 (meaning any improvement is sufficient to add the next cluster within the greedy framework)
#' @param sibling_prune numeric between 0-1 that sets the threshold for how close apparent correlation and test correlation must be for a k-cv iteration to contribute its siblings to the final chosen siblings. In other words, a lower number is more stringent, since it means the overfitting had to be really low in a k-cv iteration for it to contribute to the final sibling output.
#' @param k the k parameter in k fold cross validation (i.e. train/test partitions). Default is 5
#' @param condensed_output return output that is condensed and summarized across k_cv iterations, specifically with regard to feature importance 
#'
#' @return returns the list given by get_top_clusters for each n fold k cv run and includes the test correlation and train/test splits from each iteration
#' @export
#'
#' @examples
#' C <- 100  # represents samples
#' R <- 200 # represents features
#' y <- rnorm(C) # represents outcome variable
#' X <- matrix(rbeta(R*C, 2, 3), nrow = R)  # simulate data matrix
#' cv_res <- run_cross_validation
#'
run_cross_validation <- function( y,
                                  X,
                                  Verbose.pass = FALSE,
                                  restrict = FALSE,
                                  rank = FALSE,
                                  Verbose = TRUE,
                                  standardize_features = TRUE,
                                  cluster_corr_prop = 1,
                                  ct = 1.0,
                                  sibling_prune = 0.10,
                                  k = 5,
                                  condensed_output = TRUE
){
  
  if(isTRUE(Verbose)){
    print(paste0('Running ', k, '-fold cross validation'), quote = F)
  }
  
  # divide y into k parts
  splits <- pact::KfoldCV(length(y), k)
  k_cv_res <- sapply(1:k, function(x){
    
    if (isTRUE(Verbose)){
      print(paste0('    iteration ', x), quote = F)
    }
    
    # split data into test and train
    y.train <- y[splits != x]
    y.test <- y[splits == x]
    X.train <- X[, splits != x]
    X.test <- X[, splits == x]
    
    if (isTRUE(Verbose)){
      print('    Calculating TSCs', quote = F)
    }
    
    # run kTSCR
    TSC.res <- get_top_clusters(y = y.train,
                              X = X.train,
                              Verbose = Verbose.pass,
                              restrict = restrict,
                              rank = rank,
                              standardize_features = standardize_features,
                              cluster_corr_prop = cluster_corr_prop,
                              ct = ct
    )
    
    if (isTRUE(Verbose)){
      print('    Calculating K Scores', quote = F)
    }
    
    # use pairs from train data to get K score from test data
    K.test <- calc_K_from_pairs(TSC.res$SiblingIndices, X.test, restrict = restrict)
    
    if (isTRUE(Verbose)){
      print('    Testing Correlation', quote = F)
    }
    
    # get correlation between K.test and y.test
    test.cor <- suppressWarnings(stats::cor(y.test, K.test)) # suppress warnings that occur when there is no variation across K in the test set
    
    if (is.na(test.cor)){
      if (all_equal(K.test)){
        # adding jitter so correlation doesn't produce NA. This should still lead to a POOR cor value, which is what I want in this case.
        K.test <- K.test + stats::rnorm(length(K.test))
        test.cor <- stats::cor(y.test, K.test)
      } else{
        print("ERROR - is.na(test.cor) but !all_equal(K.test)")
        stop()
      }
    }
    
    app.cor <- TSC.res$ApparentCorrelation
    
    feature_importance <- sort(table(TSC.res$Siblings), decreasing = TRUE)

    return(list(app_cor = app.cor, test_cor = test.cor, elders = TSC.res$Elders, siblings = TSC.res$Siblings,feature_importance = feature_importance))
    
  })
  
  if (isTRUE(condensed_output)){
    return(condense_k_cv_output(X, k_cv_res, k, app_minus_test_thresh = sibling_prune))
  } else{
    return(k_cv_res)
  }
  
}

