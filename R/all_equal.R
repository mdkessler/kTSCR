#' Check vector items for equality
#'
#' Tests whether every element in a vector is equal
#'
#' @param vec a vector of any type that supports checking for equality (e.g. character, numeric, etc)
#'
#' @return logical
#' @export
#'
#' @examples
#' all_equal(rep("apple",10))
#' all_equal(rep(c("apple","pear"),10))
#'
all_equal <- function(vec){
  # Input: vector
  
  # Output: single boolean value (True of False)
  
  # Action: tests whether every element in vec is
  # equal
  
  return(all(vec == vec[1]))
}
