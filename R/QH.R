#' QH: Quesenberry and Hurst Method for Simultaneous Confidence Intervals
#'
#' Computes simultaneous confidence intervals for multinomial proportions
#' using the Quesenberry and Hurst (1964) method. The function calculates
#' lower and upper confidence limits for each category, adjusts them to remain
#' within the valid [0, 1] range, and computes the overall interval volume
#' (the product of the interval widths).
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
#' proposed by Quesenberry and Hurst (1964) for multinomial proportions.
#' It adjusts each interval to ensure limits remain within the [0, 1] range.
#'
#' @references
#' Quesenberry, C. P., and Hurst, D. C. (1964).
#' *Large Sample Simultaneous Confidence Intervals for Multinomial Proportions.*
#' Technometrics, **6**, 191–195.
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
#' QH(y, z)
#'
#' @importFrom stats qchisq
#' @export
QH <- function(inpmat, alpha) {
  if (any(inpmat < 0)) {
    stop("All elements of 'inpmat' must be non-negative integers.")
  }
  
  k <- length(inpmat)
  s <- sum(inpmat)
  
  if (s == 0) {
    stop("Sum of 'inpmat' must be positive.")
  }
  
  chi <- qchisq(1 - alpha, df = k - 1)
  pi <- inpmat / s
  
  QH.UL <- (chi + 2 * inpmat + sqrt(chi^2 + 4 * inpmat * chi * (1 - pi))) / (2 * (chi + s))
  QH.LL <- (chi + 2 * inpmat - sqrt(chi^2 + 4 * inpmat * chi * (1 - pi))) / (2 * (chi + s))
  
  LLA <- ifelse(QH.LL < 0, 0, QH.LL)
  ULA <- ifelse(QH.UL > 1, 1, QH.UL)
  
  diA <- ULA - LLA
  VOL <- round(prod(diA), 8)
  
  # Combine results into a data frame using final convention
  result <- data.frame(
    mean      = pi,
    lower_Lt  = LLA,
    upper_Lt  = ULA,
    Width     = diA,
    volume    = VOL
  )
  
  return(result)
  
  
}
