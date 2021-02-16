#' Predict y from K
#'
#' @param y vector with quantitative outcome variable
#' @param K K score calculated by kTSCR procedure
#'
#' @return numeric vector of predicted values for y
#' @export
#'
#' @examples
#' 
predict <- function(y, K){

  # TSP.res: list object returned by getTopPairs()

  # uses covariance to predict y according to the equation
  # (y - E(y)) / sd(y) = cov(k, y) * (K - E(K)) / sd(K), ==
  # E(y) = y - (cov(k, y) * (K - E(K)) * sd(y) / sd(K))
  # y = E(y) + (cor(k, y) * (K - E(K)) * sd(y) / sd(K))


  # error check - TODO...make test
  #stopifnot(length(y) == length(K))

  # y.pred <- mean(TSP.res$y) -
  #            (cov(K,y) * (K - mean(K)) * sd(y) / sd(K))

  y.pred <- mean(y) +
    (stats::cor(K, y) *
       (K - mean(K)) *
       stats::sd(y) /
       #sqrt(TSP.res$Variance))
       stats::sd(K))

  return(y.pred)
}
