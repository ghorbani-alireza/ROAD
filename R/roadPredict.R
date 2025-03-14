#' Predict Using a Fitted ROAD Model
#'
#' This function predicts class labels for new data using a fitted ROAD model.
#'
#' @param xtest A numeric matrix of test data (n x p).
#' @param fit A fitted ROAD model object.
#' @param fitCV A cross-validated ROAD model object.
#' @return A list containing:
#'   \item{class}{Predicted class labels for the optimal lambda.}
#'   \item{classAll}{Predicted class labels for all lambda values.}
#' @export
roadPredict <- function(xtest, fit, fitCV) {
  w <- fitCV$w
  mua <- fit$mua
  wPath <- fit$wPath
  
  ntest <- nrow(xtest)
  if (fit$para$sRoad != 0) {
    xtest <- xtest[, fit$para$sInd]
  }
  
  # Debugging: Print dimensions and values
  cat("Dimensions of xtest:", dim(xtest), "\n")
  cat("Dimensions of w:", length(w), "\n")
  cat("Dimensions of mua:", length(mua), "\n")
  cat("Dimensions of wPath:", dim(wPath), "\n")
  
  # Calculate predictions
  class <- (xtest - matrix(rep(mua, ntest), ntest, byrow = TRUE)) %*% w > 0
  classAll <- (xtest - matrix(rep(mua, ntest), ntest, byrow = TRUE)) %*% wPath > 0
  
  # Debugging: Print predictions
  cat("Class predictions:", head(class), "\n")
  cat("ClassAll predictions:", head(classAll), "\n")
  
  return(list(class = class, classAll = classAll))
}