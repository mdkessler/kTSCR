# kTSCR implementation

# functions

# addresses rounding error - if value is small enough, it makes it zero so that sqrt doesn't return NA
# special.sqrt <- function(x){ #FILL IN
#   ifelse(x < -1e-10, NaN, sqrt(pmax(x, 0)))
# }

AllEqual <- function(vec){
  # Input: vector

  # Output: single boolean value (True of False)

  # Action: tests whether every element in vec is
  # equal

  return(all(vec == vec[1]))
}

# MatMean <- function(mat.list){
#   # Input: list of matrices
#   
#   # Output: matrix
#   
#   # Action: calculates the mean value of each cell
#   # (e.g. Mij) across all matrices in mat.list.
#   # Returns these mean values in a single matrix
#   
#   # require that each matrix in mat.list is of the same
#   # size
#   stopifnot(AllEqual(sapply(mat.list, dim)[1,])) # require same num of rows
#   stopifnot(AllEqual(sapply(mat.list, dim)[2,])) # require same num of cols
#   
#   # make list into array
#   mat.cbind <- do.call(cbind, mat.list)
#   mat.arr <- array(mat.cbind, dim=c(dim(mat.list[[1]]), length(mat.list)))
#   
#   return(colMeans(aperm(mat.arr, c(3,1,2)), na.rm = TRUE))
# }

# getPairwiseMatrices <- function(X, restrict = FALSE){
#   # Input: feature matrix X
#   
#   # Output: list of matrices, with each matrix
#   # containing a score reflecting a pairwise
#   # comparison between features, per sample (i.e.
#   # each matrix reflects pairwise feature scores per
#   # sample)
#   
#   # Optional:
#   # restrict - Defaults to FALSE. Otherwise, must
#   # be a vector of feature names to restrict
#   # analysis to. If provided, the output matrices
#   # will only contain pairwise feature scores for
#   # comparisons that always include at least one
#   # of the features from the restrict vector. In
#   # other words, the rows of the output matrices
#   # will be restricted to the features in restrict
#   # although the cols will still include all features
#   # This way, all restrict features are compared to 
#   # all other features, which allows for the
#   # selection of non-restrict features as top pairs,
#   # which may reflect important feature relationships.
#   # Otherwise, if restrict == false, all features
#   # will be considered
#   
#   # Note: This function represents the scoring
#   # function by which feature pairs are chosen.
#   # That is, this function implements the indicator
#   # function that compares Xi < Xj across all samples
#   # and all genes
#   # MODIFY/REPLACE this function implementing a
#   # different scoring function
#   
#   # Action: this function compares the magnitude of
#   # features in pairwise fashion (e.g Xi < Xj).
#   # This is done across all genes per sample, and
#   # generates a list of num_samples matrices,
#   # where each matrix is an num_rows x num_rows matrix
#   # with a 1 whereever feature has a bigger magnitude
#   # than other features
#   
#   # first initializes zero matrix with dimensions num_rows x num_rows...used below to generate lower.triangle matrix
#   S <- matrix(0, nrow = nrow(X), ncol = nrow(X)) # init r x r matrix
#   
#   # generate lower triangle matrix
#   LT <- lower.tri(S) * 1 # lower triangle matrix - used below for rank-based feature comparison
#   
#   # rank features per sample - used to compare feature magnitudes in pairwise fashion
#   ranks <- apply(X, 2, rank) # feature ranks per sample
#   ranks <- split(ranks, rep(1:ncol(ranks), each = nrow(ranks))) # splits each col into a list
#   
#   # use LT and ranks to generate a list of pairwise comparison matrices (one per sample)
#   Is <- lapply(ranks, function(x){ # list of I matrices, where each I matrix is an R x R matrix with 1 whereever feature has a bigger magnitude than other features
#     if (isFALSE(restrict)){
#       LT[x, x]
#     } else{
#       LT[x, x][rownames(X) %in% restrict, ]
#     }
#   })
#   
#   # return
#   return(Is)
# }

getPairwiseRankMatrices <- function(X, restrict = FALSE, UseDiffRank = T){
  # Input: feature matrix X
  
  # Output: list of matrices, with each matrix
  # containing a score reflecting a pairwise
  # comparison between features, per sample (i.e.
  # each matrix reflects pairwise feature scores per
  # sample)
  
  # Optional:
  # restrict - Defaults to FALSE. Otherwise, must
  # be a vector of feature names to restrict
  # analysis to. If provided, the output matrices
  # will only contain pairwise feature scores for
  # comparisons that always include at least one
  # of the features from the restrict vector. In
  # other words, the rows of the output matrices
  # will be restricted to the features in restrict
  # although the cols will still include all features
  # This way, all restrict features are compared to 
  # all other features, which allows for the
  # selection of non-restrict features as top pairs,
  # which may reflect important feature relationships.
  # Otherwise, if restrict == false, all features
  # will be considered
  
  # Note: This function represents the scoring
  # function by which feature pairs are chosen.
  # That is, this function saves the rank difference between
  # Xi < Xj across all samples
  # and all features
  # MODIFY/REPLACE this function implementing a
  # different scoring function
  
  # Action: this function compares the magnitude of the rank of the
  # features in pairwise fashion (e.g Xi < Xj).
  # This is done across all genes per sample, and
  # generates a list of num_samples matrices,
  # where each matrix is an num_rows x num_rows matrix
  # with a rank difference between features
  
  # first initializes zero matrix with dimensions num_rows x num_rows...used below to generate lower.triangle matrix
  S <- matrix(0, nrow = nrow(X), ncol = nrow(X)) # init r x r matrix
  
  # generate lower triangle matrix
  LT <- lower.tri(S) * 1 # lower triangle matrix - used below for rank-based feature comparison
  
  # rank features per sample - used to compare feature magnitudes in pairwise fashion
  ranks <- apply(X, 2, rank) # feature ranks per sample
  ranks <- split(ranks, rep(1:ncol(ranks), each = nrow(ranks))) # splits each col into a list
  
  if (isTRUE(UseDiffRank)){
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


# getEIyt <- function(Is, yt){
#   # Input: 1) list of matrices, each with pairwise
#   # feature scores 2) centered outcome variable
#   
#   
#   
#   # Output: matrix
#   
#   # Action: returns a matrix with values representing the expectation (E)
#   # (across samples) of pairwise comparisons (I) x y tilde
#   # (E[I * yt])
#   # yt is the centered quantitative variable which we are
#   # looking to predict
#   
#   Iyts <- sapply(1:length(Is), function(x){  # X (i.e. I(Xi<Xj)) * yt
#     Is[[x]]*yt[x]
#   }, simplify = F)
#   
#   # To get expectation, take the mean of these matrices across samples
#   # note for clairty: XQs is already a list of matrices - one per sample - MatMean facilitates taking the mean matrices across these matrices
#   EIyt <- MatMean(Iyts)
#   
#   # return
#   return(EIyt)
# }
# 
# getPij <- function(Is, C, pseudocounts = T){
#   
#   # calculate the expectation of Xi<Xj for all pairwise comparisons
#   # this expectation can be estimated as a proportion across samples
#   # C: num cols
#   # pseudocounts: if T (default), adds small values to
#   # numerator and denominator to address divide by
#   # zero errors when Pij*(1-Pij) == 0. Otherwise, if
#   # pseudocounts set to FALSE, divide by zero issue is
#   # addressed by setting all Pij*(1-Pij) == 0 <- NA
#   
#   if (isTRUE(pseudocounts)){
#     Pij <- (Reduce('+', Is)  + 1) / (C + 2) # proportion Xi<Xj across samples for all comparisons
#   } else{
#     Pij <- Reduce('+', Is) / C # proportion Xi<Xj across samples for all comparisons
#     Pij[Pij*(1-Pij) == 0] <- NA
#   }
#   
#   # return
#   return(Pij)
# }
# 
# getEIKt <- function(Is, Kt){
#   # calculate Is * Kt
#   # this represents a multiplication between the
#   # pairwise comparison matrices and the centered K score value,
#   # per sample
#   
#   # list of matrices representing I(Xi<Xj) * Kt, per sample
#   IKs <- sapply(1:length(Is), function(x){  # I(Xi<Xj) x Ktilde
#     Is[[x]]*Kt[x]
#   }, simplify = F)
#   
#   # get the mean matrix across XKs matrices
#   EIKt <- MatMean(IKs)
#   
#   # return
#   return(EIKt)
# }

MakeIndicatorFeatures <- function(Is){
  num_features <- nrow(Is[[1]]) * (nrow(Is[[1]]) - 1) / 2
  #num_features <- nrow(Is[[1]])**2 - nrow(Is[[1]])
  num_samples <- length(Is)
  I_mat <- matrix(rep(0, num_samples*num_features),
                  nrow = num_features)
  
  new_var_names <- rep(NA, (nrow(Is[[1]])-1 * nrow(Is[[1]]) / 2))
  #new_var_names <- rep(NA, nrow(Is[[1]])**2 - nrow(Is[[1]]))
  feature_counter <- 1
  for (i in 1:(nrow(Is[[1]])-1)){
  #for (i in 1:nrow(Is[[1]])){
    for (j in (i+1):nrow(Is[[1]])){
    #for (j in 1:nrow(Is[[1]])){
      # if (i == j){
      #   next
      # }
      for (l in 1:length(Is)){
        I_mat[feature_counter,l] <- Is[[l]][j,i]
      }
      # set new var name
      #new_var_name <- paste0('V', j, '_', 'V', i) # Delete me when refactoring
      new_var_name <- paste0(j, '_', i) # Delete me when refactoring
      # save new var name
      new_var_names[feature_counter] <- new_var_name
      # increment counter
      feature_counter <- feature_counter + 1
    }
  }
  
  # transpose matrix which sets variables as cols
  I_mat <- t(I_mat)
  # set column names of I_mat with new_var_names
  colnames(I_mat) <- new_var_names
  
  # if any column is entirely identical at every element,
  # add jitter to prevent issues when calculating correlations with y
  I_mat <- apply(I_mat, 2, function(col){
    if(AllEqual(col)){
      return(col + rnorm(length(col)))
    }else{
      return(col)
    }
  })
  
  return(I_mat)
}

SortedCorrIFeatures <- function(I_mat, y){
  # Return a sorted list of indicator features based on their pearson correlation with y 
  sorted_corrs <- sort(apply(I_mat, 2, function(x){
    cor(x,y)
    }), decreasing = T)
}

# NOTE (FIX): Gets Elders and Brothers - refactor this
GetBrothers <- function(sorted_corrs, elder_index, cc = 3){
  # choose the pairs that include the top feature hit
  # these other features in these chosen pairs are called brothers
  # and they define the "cluster" for this top feature
  # cc: cluster contour - a parameter for how many consecutive misses until I stop saving cluster brothers
  # FILL IN addition description
  
  var_list <- list()
  for (var in unique(unlist(sapply(names(sorted_corrs), strsplit, "_")))){
    var_list[var] <- 0
  }
  
  for (i in 1:length(sorted_corrs)){
    var1_var2 <- names(sorted_corrs)[i]
    corr <- abs(as.numeric(sorted_corrs[i]))
    var1 <- str_split_fixed(var1_var2, "_", n = 2)[1]
    var2 <- str_split_fixed(var1_var2, "_", n = 2)[2]
    var_list[[var1]] <- var_list[[var1]] + min(abs(1/log(corr)),100)
    var_list[[var2]] <- var_list[[var2]] + min(abs(1/log(corr)),100)
  }
  
  # the top hit is called the elder
  sorted_var_list <- sort(unlist(var_list), decreasing = T)
  elder_val <- sorted_var_list[[elder_index]]
  elder <- names(sorted_var_list)[elder_index]
  
  #last_hit <- 0 # counts consecutive misses without a top hits - i.e. the number of previous top Indicator features without a top_hit
  elder_val_build_up <- 0 # init at zero - will record a sum of how much of the elder value has been found as every brother pair is found
  vars_to_keep <- c() # stores brother vars to keep
  for (i in 1:length(sorted_corrs)){
    var1_var2 <- names(sorted_corrs)[i]
    var1 <- str_split_fixed(var1_var2, "_", n = 2)[1]
    var2 <- str_split_fixed(var1_var2, "_", n = 2)[2]
    # reverse the Var1_Var2 order if the correlation is negative
    corr <- as.numeric(sorted_corrs[i])
    if (corr < 0){
      var1_var2 <- paste(var2, var1, sep = "_")
    }
    
    # break if you missed the allowed number of consecutive pairs without an elder
    #if (last_hit == cc){ # cluster contour - controls how many brothers are included based on how many consecutive skips are allowed
    #  break
    if (elder_val_build_up / elder_val >= cc){ # cluster contour - controls how many brothers are included based on how many consecutive skips are allowed
      break
    
    } else{ # otherwise, see if current pair has an elder in it - if so, add it to vars to keep, which will define brothers
      if (var1 == elder | var2 == elder){
        vars_to_keep <- append(vars_to_keep, var1_var2)
        #last_hit <- 0 # reset consecutive top hits counter
        # add cor value to elder_val_build_up
        elder_val_build_up <- elder_val_build_up + abs(corr)
      } # else{
      #   last_hit <- last_hit + 1 # no top hit
      # }
    }
  }
  
  return(list(elder, vars_to_keep))
  
}

GetBrothersIndices <- function(brothers){
  # convert brother pair feature names into indices from the
  # original feature matrix X
  
  indices <- t(apply(do.call(rbind, strsplit(gsub("V", "", brothers), "_")), 1, as.numeric))
  
  return(indices)
}


CalcKfromPairs <- function(pairs, X, restrict){
  
  # calculate K score given feature matrix and a set of
  # identified pairs
  # X is feature matrix to use for calculating K Score
  
  # convert X to Is
  Is <- getPairwiseRankMatrices(X, restrict)
  
  KScore <- unlist(lapply(Is, function(j){ # iterate over matrices in Is list
    # j is a matrix from Is (represents a samples rank diff matrix)
    diffs <- apply(pairs, 1, function(i){ # iterates over feature pair indices
      return(j[i[1], i[2]]) # the way the pairs
    })
    return(sum(diffs))
  }))
  
  return(KScore)
  
}

ReduceIMat <- function(I_mat, elder){
  # Remove all columns from I_mat that have elder as one of the
  # indicator comparison pairs
  
  exclude <- which(sapply(1:ncol(I_mat), function(i){
    elder %in% unlist(strsplit(colnames(I_mat)[i], "_"))
  }))
  
  # deal with the special case where there is only one indicator feature column left in the reduced (i.e. updated) matrix
  if (ncol(I_mat) - length(exclude) == 1){
    ret <- as.matrix(I_mat[, -exclude]) # as matrix is necessary to explictly convert 1 column matrix from vector back to matrix
    colnames(ret) <- colnames(I_mat)[!(1:ncol(I_mat)) %in% exclude]
    return(ret)
  }
  # else
  return(I_mat[, -exclude]) 
}

getTopClusters <- function(y,
                        X,
                        #kmax = nrow(X) / 2,
                        Verbose = F,
                        restrict = F,
                        rank = F,
                        standardize_features = T,
                        cc = 1,
                        ct = 1.0, # 1.01 would mean that correlation must improve by more than 1% to keep searching for elders/clusters during greedy optimization
                        UseDiffRank = T
                        ){
  
  # y: quantitative data
  # X: feature matrix, samples as cols
  # kmax: max value of top pairs to output
  # if kmax is not set, set it to the number of
  # features divided by 2 (as this is the
  # theoretical max if every feature was in a pair
  # once and only once)
  # rank: if rank is T, rank(y) is used instead of y
  # Verbose: whether to print verbosely
  # restrict: features to restrict to
  # standardize_features: boolean indicating whether
  # features should be standardized
  
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
  
  # PERFORM INITIAL CALCULATIONS
  
  # center y
  if (Verbose){
    print('Centering outcome values (i.e. y)', quote = F)
  }
  # yt <- scale(y, scale = F, center = T) # y tilde - centers y
  
  # generate list of  pairwise comparison matrices
  if (Verbose){
    print('Calculating pairwise indicator matrices', quote = F)
  }
  
  Is <- getPairwiseRankMatrices(X, restrict, UseDiffRank)
  I_mat <- MakeIndicatorFeatures(Is)
  
  # # calculate E[X * yt]
  # if (Verbose){
  #   print('Calculating E[X * yt] matrix', quote = F)
  # }
  # EIyt <- getEIyt(Is, yt)
  # 
  # # calculate proportion I(Xi<Xj) across all samples
  # if (Verbose){
  #   print('Calculating Pij matrix', quote = F)
  # }
  # Pij <- getPij(Is, ncol(X)) 
  # 
  # # INITIALIZE VARIABLES
  # 
  if (Verbose){
  
    print('Initializing variables', quote = F)
  }

  K <- rep(0, ncol(X))  # init K <- 0 for every sample (stores running K score)
  G <- list() # init list of top hit features (termed Elders)
  B <- list() # init list to hold list of brothers for each eldest
  k <- 1 # little k - pair you're up to - max K must be <= k each iteration
  last_corr <- 0
  # use a while loop to identify pairs that optimize correlation

  if (Verbose){
    print('Starting kTSCR iteration', quote = F)
  }
  
  while(TRUE){
    # if kmax is provided as argument, end loop if k > kmax
    # remember, the way k is set, it is always +1 past the number
    # of pairs already stored
    # if (k > kmax){ # stop after k pairs
    #   k <- k - 1 # decrease k to match num pairs chosen
    #   break # break out of while loop
    # }

    if (Verbose){
      print(paste0("    Iteration ", k), quote = F)
    }

    # # the following code makes sure I don't include
    # # the same features in more than one pair
    # # In other words, if a feature was already chosen
    # # as part of a previous pair, then the following
    # # code sets the values of this feature to NA so
    # # that it can't be chosen by optimization procedure
    # # if (length(G) > 0){ # if G already has pairs, (i.e. beyond the first loop iteration)
    # #   EIyt[last_pair[1], ] <- NA # set row from row index in last pair to NA
    # #   EIyt[, last_pair] <- NA # set cols matching indices in last pair to NA
    # #   if (last_pair[2] <= nrow(EIyt)){ # set row of second pair to NA if the index exists, which it may not given that you may have restriced the feature comparisons
    # #     EIyt[last_pair[2], ] <- NA # set row from col index in last pair to NA
    # #   }
    # # }
    # # allow features to repeat across pairs
    # # if (length(G) > 0){ # if G already has pairs, (i.e. beyond the first loop iteration)
    # #   EIyt[last_pair[1], last_pair[2]] <- NA # set row from row index in last pair to NA
    # #   if (last_pair[2] <= nrow(EIyt)){ # set row of second pair to NA if the index exists, which it may not given that you may have restriced the feature comparisons
    # #     EIyt[last_pair[2], last_pair[1]] <- NA # set row from col index in last pair to NA
    # #   }
    # # }
    # #
    # 
    # # scale updated K score
    # Kt <- scale(K, scale = F, center = T) # K Tilde
    # 
    # EIKt <- getEIKt(Is, Kt) # update E[X *Kt]
    # 
    # # Again, makes sure I don't include the same features in more than one pair
    # # if (length(G) > 0){ # if G already has pairs, (i.e. beyond the first loop iteration)
    # #   EIKt[last_pair[1], ] <- NA
    # #   EIKt[, last_pair] <- NA
    # #   Pij[last_pair[1], ] <- NA
    # #   Pij[, last_pair] <- NA
    # #   if (last_pair[2] <= nrow(EIyt)){ # set row of second pair to NA if the index exists, which it may not given that you may have restriced the gene comparisons
    # #     EIKt[last_pair[2], ] <- NA
    # #     Pij[last_pair[2], ] <- NA
    # #   }
    # # }
    # # allow features to repeat across pairs
    # # if (length(G) > 0){ # if G already has pairs, (i.e. beyond the first loop iteration)
    # #   EIKt[last_pair[1], last_pair[2]] <- NA
    # #   Pij[last_pair[1], last_pair[2]] <- NA
    # #   if (last_pair[2] <= nrow(EIyt)){ # set row of second pair to NA if the index exists, which it may not given that you may have restriced the gene comparisons
    # #     EIKt[last_pair[2], last_pair[1]] <- NA
    # #     Pij[last_pair[2], last_pair[1]] <- NA
    # #   }
    # # }

    # # calculate the term you will maximize
    # numer <- c + EIyt # numerator of term
    # denom <- special.sqrt(v + (2 * EIKt) + (Pij*(1-Pij))) # denominator of term
    # arg_max_mat <- numer/denom # values of term to maximize in matrix
    # 
    # # set diagonal of arg_max_mat to NA so no pairs with the same
    # # genes can be chosen
    # diag(arg_max_mat) <- NA

    # # choose feature pair that maximizes numer/denom
    # max_corr_val <- max(arg_max_mat, na.rm = T) # max value
    # # If no pairs can increase correlation, break pair selection
    # if ((max_corr_val/sqrt(var(y))) < last_corr){ # must normalize max_corr_cval by sqrt(var(y)) to put it in units of proper correlation
    #   if (Verbose){
    #     print("    No pairs can increase correlation - BREAKING ITERATION")
    #   }
    #   break
    # }
    # 
    # # Find cells with max values
    # max_cells <- which(arg_max_mat == max_corr_val, arr.ind = TRUE) # indexes of all cells that match the max value
    # 
    # # FILL IN - HANDLE TIES HERE - UPDATE: see below, where random rop pair choice is implemented to handle ties
    # # For now, just print a warning!
    # # if (dim(max_cells)[1] != 1 | dim(max_cells)[2] != 2){
    # #   #warning(dim(max_cells))
    # #   print(paste("ERROR - dim(max_cells):", dim(max_cells)))
    # #   print(paste("max_cells",max_cells))
    # #   print(paste("max_corr_val",max_corr_val))
    # #   #View(arg_max_mat)
    # #   print(paste("c/sqrt(v)",c/sqrt(v)))
    # #   print(paste("EIyt[1,1]",EIyt[1,1]))
    # #   print(paste("2 * EIKt[1,1]",2 * EIKt[1,1]))
    # #   print(paste("Pij[1,1]*(1-Pij[1,1])",Pij[1,1]*(1-Pij[1,1])))
    # #   stop()
    # # }
    # 
    # # break ties randomly
    # 
    # pair_choice_idx <- sample(1:dim(max_cells)[1], 1) # if there are ties, this will choose one pair randomly. If not, this will choose the top pair
    # i <- max_cells[pair_choice_idx,1] # row corresponding to max value
    # j <- max_cells[pair_choice_idx,2] # col corresponding to max value
    
    # if all features have already been found to contain a brother, break
    if (dim(I_mat)[2] == 0){
      if (Verbose){
        print(paste0("Converged at k=", k,": ALL BROTHERS"), quote = F)
      }
      break
    }
    
    # get correlations
    if (k == 1){ # first iteration
      sorted_corrs <- SortedCorrIFeatures(I_mat, y)
    } else{
      y_err <- y# - K
      sorted_corrs <- SortedCorrIFeatures(I_mat, y_err)
      if (dim(I_mat)[2] == 1){ #FIX
      }
    }
    
    # find elders and brothers
    brothers <- c()
    elder_index <- 0 # increments with each pass of while loop
    while (length(brothers) == 0){
      elder_index <- elder_index + 1
      get_brothers <- GetBrothers(sorted_corrs, elder_index, cc = cc)
      # extract elder
      elder <- get_brothers[[1]]
      # extract brothers 
      brothers <- get_brothers[[2]]
      
      # if this isn't the first pair, then we know we at least have some brothers from the previous iteration
      if (k > 1){
        break
      }
      
    }
    
    # Test for a WHILE LOOP BREAK CONDITION -
    # if no brother are found, break the while loop and DO NOT
    # save the elder or any updates to any variables
    # i.e. the greedy optimization has finished
    if (is.null(brothers)){
      if (Verbose){
        print(paste0("Converged at k=", k,": NO BROTHERS"), quote = F)
      }
      break
    }
    
    # if you don't break...
    # convert brother pairs into indices for use in calculating K score with them
    brothers <- GetBrothersIndices(brothers)
    # Test for another WHILE LOOP BREAK CONDITION - CORRELATION CONVERGENCE
    new_K <- K + CalcKfromPairs(brothers, X, restrict = restrict)
    if (sd(new_K) == 0 | sd(y) == 0){
      print(new_K)
      print(y)
    }
    current_corr <- cor(new_K, y)
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
    B[[k]] <- brothers
    
    # update k for next iteration (remember, in this implementation, k is one ahead whenever I break the while loop)
    k <- k + 1
    
    # update last_corr
    last_corr <- current_corr
    
    # update I_mat
    I_mat <- ReduceIMat(I_mat, elder)
  }
  
  # OUTPUT RELEVANT INFORMATION

  # Now that you are outside the while loop,
  # return relevant information as a list
  # Note: look further into customized object to return data in

  # store indices of elders and brother
  elder_indices <- as.numeric(unlist(G))
  brother_indices <- do.call(rbind, B)
  
  if (!is.null(rownames(X))){ #& !ReturnPairsAsIndices){
    elders <- sapply(G, function(v){ # return vector of elders
                        # i <- as.numeric(gsub(pattern = "V", replacement = "", x = v)) # delete me when refactoring
                        i <- as.numeric(v)
                        rownames(X)[i]
                      })
    brothers <- do.call(rbind, # return matrix of pairs (each row is a brother pair)
                        lapply(B, function(v){
                          # i <- as.numeric(gsub(pattern = "V", replacement = "", x = v)) # delete me when refactoring
                          i <- as.numeric(v)
                          matrix(rownames(X)[i], byrow = F, nrow = nrow(v))
                        }))
  } else{
    elders <- NULL
    brothers <- NULL
  }
  
  # list to output
  outlist <- list(y = y,
                  KScore = K,
                  Correlation = cor(K, y),
                  ElderIndices = elder_indices,
                  BrotherIndices = brother_indices,
                  Elders = elders,
                  Brothers = brothers
                  )
  
  # DELETE - init_arg_max = init_arg_max)

  # return
  return(outlist)

  if (Verbose){
    print('DONE', quote = F)
  }
}

# kTSCpredict <- function(TSP.res){
#   
#   # TSP.res: list object returned by getTopPairs()
#   
#   # uses covariance to predict y according to the equation
#   # (y - E(y)) / sd(y) = cov(k, y) * (K - E(K)) / sd(K), ==
#   # E(y) = y - (cov(k, y) * (K - E(K)) * sd(y) / sd(K))
#   # y = E(y) + (cor(k, y) * (K - E(K)) * sd(y) / sd(K))
#   
#   
#   # error check
#   stopifnot(length(TSP.res$y) == length(TSP.res$KScore))
#   
#   # y.pred <- mean(TSP.res$y) - 
#   #            (cov(K,y) * (K - mean(K)) * sd(y) / sd(K))
#   
#   y.pred <- mean(TSP.res$y) +
#     #(TSP.res$Covariance *
#     (cor(TSP.res$KScore, TSP.res$y) *
#        (TSP.res$KScore - mean(TSP.res$KScore)) *
#        sd(TSP.res$y) /
#        sqrt(TSP.res$Variance))
#   
#   
#   return(y.pred)
# }

runCV <- function(y,
                  X,
                  Verbose.pass = F,
                  restrict = F,
                  rank = F,
                  Verbose = T,
                  standardize_features = T,
                  cc = 1,
                  ct = 1.0,
                  k_cv = 5
                  
){
  
  # k_cv : k fold for cross validation
  # other options are to be passed through to getTopPairs()
  # note: Verbose.pass gets passed to getTopPairs()
  # whereas Verbose controls verbosity of this runCV function
  
  if(isTRUE(Verbose)){
    print(paste0('Running ', k_cv, '-fold cross validation'), quote = F)
  }
  
  # divide y into k parts
  require(pact)
  splits <- pact::KfoldCV(length(y), k_cv)
  k_cv_cors <- sapply(1:k_cv, function(x){
    
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
    
    # if correlation filtering is true, add it here!
    # if you implement it, add verbosity message
    
    # run kTSCR
    TSC.res <- getTopClusters(y = y.train,
                           X = X.train,
                           Verbose = Verbose.pass,
                           restrict = restrict,
                           rank = rank,
                           standardize_features = standardize_features,
                           cc = cc,
                           ct = ct
                           )
    
    if (isTRUE(Verbose)){
      print('    Calculating K Scores', quote = F)
    }
    
    # use pairs from train data to get K score from test data
    K.test <- CalcKfromPairs(TSC.res$BrotherIndices, X.test, restrict = restrict)
    
    if (isTRUE(Verbose)){
      print('    Testing Correlation', quote = F)
    }
    
    # get correlation between K.test and y.test
    test.cor <- suppressWarnings(cor(y.test, K.test)) # suppress warnings that occur when there is no variation across K in the test set

    if (is.na(test.cor)){
      if (AllEqual(K.test)){
        # adding jitter so correlation doesn't produce NA. This should still lead to a POOR cor value, which is what I want in this case.
        K.test <- K.test + rnorm(length(K.test))
        test.cor <- cor(y.test, K.test)
      } else{
        print("ERROR - is.na(test.cor) but !AllEqual(K.test)")
        stop()
      }
    }
    
    # return 
    TSC.res$test.cor <- test.cor
    TSC.res$splits <- splits
    
    return(TSC.res)
    
  })
  
  return(k_cv_cors)
  
}
