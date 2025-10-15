#' FS: Fitzpatrick and Scott Method for Simultaneous Confidence Intervals
#'
#' Computes simultaneous confidence intervals for multinomial proportions
#' using the Fitzpatrick and Scott (FS) method. The function estimates the
#' lower and upper confidence limits for each category, adjusts them to
#' remain within the [0, 1] range, and calculates the overall volume
#' (product of interval widths).
#'
#' @param inpmat Integer vector of observed counts (non-negative values).
#' @param alpha Desired statistical Significance level.
#'
#' @return Prints the original and adjusted confidence intervals for each
#' category, as well as the overall interval volume. 
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 0.05
#' FS(y, z)
#'
#' @importFrom stats qnorm
#' @export
FS <- function(inpmat, alpha) {
  # Input checks
  if (any(inpmat < 0)) {
    stop("All elements of 'inpmat' must be non-negative integers.")
  }
  
  k <- length(inpmat)
  s <- sum(inpmat)
  if (s == 0) {
    stop("Sum of input counts ('inpmat') must be positive.")
  }
  
  # Z-value for confidence interval
  zval <- abs(qnorm(1 - (alpha / 2)))
  pi    <- inpmat / s
  
  # Compute Wald-type intervals
  FS.LL <- pi - (zval / (2 * sqrt(s)))
  FS.UL <- pi + (zval / (2 * sqrt(s)))
  
  # Ensure bounds are within [0,1]
  LLA <- pmax(FS.LL, 0)
  ULA <- pmin(FS.UL, 1)
  
  # Interval and volume
  diA <- ULA - LLA
  volume <- round(prod(diA), 8)
  
  # Combine results into a data frame
  result <- data.frame(
    mean      = pi,
    lower_Lt  = LLA,
    upper_Lt  = ULA,
    Width  = diA,
    volume    = volume
  )
  
  return(result)
}
