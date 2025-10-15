#' WALDCC: Wald Method with Continuity Correction for Simultaneous Confidence Intervals
#'
#' Computes Wald-type simultaneous confidence intervals for multinomial
#' proportions, incorporating a continuity correction. These intervals are
#' symmetric about the sample proportions and apply a small correction to
#' improve coverage accuracy, particularly for small samples.
#'
#' @param inpmat Integer vector of observed cell counts corresponding to a
#' categorical dataset. All values must be non-negative.
#' @param alpha Desired statistical significance level
#'
#' @return
#' Prints the original and adjusted confidence intervals for each category,
#' along with the volume (product of interval widths). 
#'
#' @details
#' The correction term 1/2n ensures more accurate interval bounds,
#' especially when the proportions are near 0 or 1.
#'
#' @references
#' Wald, A. (1943).
#' *Tests of Statistical Hypotheses Concerning Several Parameters When the Number of Observations is Large.*
#' Transactions of the American Mathematical Society, **54**, 426–482.
#'
#' @author
#' Dr. M. Subbiah
#'
#' @seealso
#' \code{\link{BMDE}}, \code{\link{WALD}}, \code{\link{SG}}
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 0.05
#' WALDCC(y, z)
#'
#' @importFrom stats qchisq
#' @export
WALDCC <- function(inpmat, alpha) {
  if (any(inpmat < 0)) stop("All elements of 'inpmat' must be non-negative.")
  if (sum(inpmat) == 0) stop("Sum of 'inpmat' must be positive.")
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be between 0 and 1.")
  
  k <- length(inpmat)
  s <- sum(inpmat)
  chi <- qchisq(1 - alpha, df = 1)
  pi <- inpmat / s
  
  WALDCC.LL <- pi - sqrt(chi * pi * (1 - pi) / s) - (1 / (2 * s))
  WALDCC.UL <- pi + sqrt(chi * pi * (1 - pi) / s) + (1 / (2 * s))
  
  LLA <- ifelse(WALDCC.LL < 0, 0, WALDCC.LL)
  ULA <- ifelse(WALDCC.UL > 1, 1, WALDCC.UL)
  
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
