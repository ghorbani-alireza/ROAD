#' Generate Logarithmically Spaced Values
#'
#' This function generates logarithmically spaced values between `10^a` and `10^b`.
#'
#' @param a Numeric scalar: the starting exponent (log10 of the first value).
#' @param b Numeric scalar: the ending exponent (log10 of the last value).
#' @param n Integer: the number of points to generate (default: 50).
#' @return A numeric vector of logarithmically spaced values.
#' @export
logspace <- function(a, b, n = 50) {
  10^seq(a, b, length.out = n)
}