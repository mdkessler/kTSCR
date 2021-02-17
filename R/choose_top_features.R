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
  km <- stats::kmeans(vec,3)
  top_features <- names(km$cluster)[km$cluster %in% as.numeric(names(sort(table(km$cluster))[1:2]))]
  
  return(top_features)
}