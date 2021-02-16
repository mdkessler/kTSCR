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
#' @param k the k parameter in k fold cross validation (i.e. train/test partitions). Default is 5
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
                                  Verbose.pass = F,
                                  restrict = F,
                                  rank = F,
                                  Verbose = T,
                                  standardize_features = T,
                                  cluster_corr_prop = 1,
                                  ct = 1.0,
                                  k = 5
){
  
  if(isTRUE(Verbose)){
    print(paste0('Running ', k, '-fold cross validation'), quote = F)
  }
  
  # divide y into k parts
  splits <- pact::KfoldCV(length(y), k)
  k_cv_cors <- sapply(1:k, function(x){
    
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
    
    feature_importance <- sort(table(c(unlist(kTSCR.cv['Elders',1:5]), unlist(kTSCR.cv['Siblings',1:5]))), decreasing = TRUE)

    return(list(app_cor = app.cor, test_cor = test.cor, feature_importance = feature_importance))
    
  })
  
  return(k_cv_cors)
  
}

