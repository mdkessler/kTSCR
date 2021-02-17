#' Choose top features
#'
#' Use kmeans clustering to separate features into priority groups based on feature importance scores calculated during n fold k-cv iterations
#'
#' @param vec a numerical vector of features importance scores. Vector should be named with feature names
#'
#' @return a vector of top_features
#' @export
#'
#' @examples
#' 
#'
choose_top_features <- function(vec){
  
  #TODO - add more regarding how k is chosen and what makes it into top features
  
  # kmeans to separate groups of prioritized features
  # you will do these with a trycatch block to handle cases where there are not enough distinct values to have three clusters
  
  km <- tryCatch({
    stats::kmeans(vec,3)
  }, warning = function(cond) {
    message(paste0("kmeans warning k = : ",3))
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
    return(NULL)
  }, error = function(cond) {
    # rerun kmeans with k = 2
      stats::kmeans(vec,2)
  }, finally=NULL)
  
  
  if (length(km$centers) == 3){
    top_features <- names(km$cluster)[km$cluster %in% as.numeric(names(sort(table(km$cluster))[1:2]))]
  } else{
    top_features <- names(km$cluster)[km$cluster %in% as.numeric(names(sort(table(km$cluster))[1]))]
  }
  
  return(top_features)
}