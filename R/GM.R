#' GM: Goodman Method for Simultaneous Confidence Intervals
#'
#' Computes simultaneous confidence intervals for multinomial proportions
#' using the Goodman (1965) method. The function calculates lower and upper
#' confidence limits for each category, adjusts them to remain within the
#' [0, 1] range, and computes the overall interval volume (the product of
#' interval widths).
#'
#' @param inpmat Integer vector of observed cell counts corresponding to
#' a categorical dataset. Must contain non-negative values.
#' @param alpha Desired statistical significance level 
#'
#' @return Prints the original and adjusted confidence intervals for each
#' category, as well as the overall volume of the simultaneous confidence
#' intervals. 
#'
#' @details
#' This function implements the simultaneous confidence interval method
#' proposed by Goodman (1965) for multinomial proportions. It adjusts each
#' interval to ensure the limits fall within the valid probability range.
#'
#' @references
#' Goodman, L. A. (1965). *On Simultaneous Confidence Intervals for Multinomial
#' Proportions.* Technometrics, **7**, 247–254.
#'
#' @author
#' Dr. M. Subbiah
#'
#' @seealso
#' \code{\link{BMDE}}, \code{\link{WALD}}, \code{\link{WS}}
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 0.05
#' GM(y, z)
#'
#' @importFrom stats qchisq
#' @export
GM <- function(inpmat, alpha) {
  # Input checks
  if (any(inpmat < 0)) {
    stop("All elements of 'inpmat' must be non-negative integers.")
  }
  
  k <- length(inpmat)
  s <- sum(inpmat)
  if (s == 0) {
    stop("Sum of 'inpmat' must be positive.")
  }
  
  # Compute chi-square threshold
  chi <- qchisq(1 - (alpha / k), df = 1)
  pi  <- inpmat / s
  
  # Goodman intervals
  GM.UL <- (chi + 2 * inpmat + sqrt(chi^2 + 4 * inpmat * chi * (1 - pi))) / (2 * (chi + s))
  GM.LL <- (chi + 2 * inpmat - sqrt(chi^2 + 4 * inpmat * chi * (1 - pi))) / (2 * (chi + s))
  
  # Ensure bounds are within [0,1]
  LLA <- pmax(GM.LL, 0)
  ULA <- pmin(GM.UL, 1)
  
  # Interval width and volume
  diA    <- ULA - LLA
  volume <- round(prod(diA), 8)
  
  # Combine results into a data frame using final convention
  result <- data.frame(
    mean      = pi,
    lower_Lt  = LLA,
    upper_Lt  = ULA,
    Width     = diA,
    volume    = volume
  )
  
  return(result)
}
