#' Calculate Root Mean Squared Error
#'
#' @param y vector of numeric values, such as the outcome variable in regression
#' @param y_pred another numeric vector, such as y values predicted by a model
#'
#' @return numeric respresent the rmse
#' @export
#'
#' @examples
#' 
rmse = function(y, y_pred){
  sqrt(mean((y - y_pred)^2))
}