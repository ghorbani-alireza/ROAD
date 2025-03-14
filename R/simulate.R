#' Simulate Data for ROAD Demonstration
#'
#' This function generates simulated data for demonstrating the ROAD algorithm. The data is generated
#' following the example from the MATLAB implementation of ROAD.
#'
#' @param p Number of variables (default: 1000).
#' @param n Number of observations per class (default: 300).
#' @param s0 Number of nonzero mean differences (default: 10).
#' @param rho Pairwise correlation coefficient (default: 0.5).
#' @param randSeed Random seed for reproducibility (default: 1).
#' @return A list containing:
#'   \item{x}{A numeric matrix of training data (n x p).}
#'   \item{y}{A numeric vector of training labels (0/1).}
#'   \item{xtest}{A numeric matrix of test data (n x p).}
#'   \item{ytest}{A numeric vector of test labels (0/1).}
#' @export
#' @examples
#' # Generate simulated data
#' sim_data <- simulate_road_data()
#' x <- sim_data$x
#' y <- sim_data$y
#' xtest <- sim_data$xtest
#' ytest <- sim_data$ytest
simulate_road_data <- function(p = 1000, n = 300, s0 = 10, rho = 0.5, randSeed = 1) {
  # Set random seed
  set.seed(randSeed)
  
  # Generate means for the two classes
  mu1 <- rep(0, p)
  mu2 <- rep(0, p)
  mu2[1:s0] <- 1
  
  # Generate the covariance matrix
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  
  # Generate training data
  Y1Train <- MASS::mvrnorm(n, mu1, Sigma)
  Y2Train <- MASS::mvrnorm(n, mu2, Sigma)
  
  # Generate testing data
  Y1Test <- MASS::mvrnorm(n, mu1, Sigma)
  Y2Test <- MASS::mvrnorm(n, mu2, Sigma)
  
  # Combine data
  x <- rbind(Y1Train, Y2Train)
  y <- c(rep(0, n), rep(1, n))
  xtest <- rbind(Y1Test, Y2Test)
  ytest <- c(rep(0, n), rep(1, n))
  
  # Return the simulated data
  list(x = x, y = y, xtest = xtest, ytest = ytest)
}