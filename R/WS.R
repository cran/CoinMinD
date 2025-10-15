#' WS: Wilson Score Method for Simultaneous Confidence Intervals
#'
#' Computes Wilson score-type simultaneous confidence intervals for multinomial
#' proportions. The Wilson method improves upon the Wald interval by ensuring
#' better coverage probabilities, especially for small samples or proportions
#' near 0 or 1.
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
#' This approach adjusts both the center and width of the confidence interval
#' to account for the sampling distribution of proportions, leading to
#' non-symmetric intervals that perform better than the simple Wald intervals.
#'
#' @references
#' Wilson, E. B. (1927).
#' *Probable Inference, the Law of Succession, and Statistical Inference.*
#' Journal of the American Statistical Association, **22**, 209–212.
#'
#' @author
#' Dr. M. Subbiah
#'
#' @seealso
#' \code{\link{BMDE}}, \code{\link{WALD}}, \code{\link{WALDCC}}, \code{\link{SG}}
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 0.05
#' WS(y, z)
#'
#' @importFrom stats qchisq
#' @export
WS <- function(inpmat, alpha) {
  if (any(inpmat < 0)) stop("All elements of 'inpmat' must be non-negative.")
  if (sum(inpmat) == 0) stop("Sum of 'inpmat' must be positive.")
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be between 0 and 1.")
  
  k <- length(inpmat)
  s <- sum(inpmat)
  chi <- qchisq(1 - alpha, df = 1)
  pi <- inpmat / s
  
  WS.UL <- (chi + 2 * inpmat + sqrt(chi^2 + 4 * inpmat * chi * (1 - pi))) / (2 * (chi + s))
  WS.LL <- (chi + 2 * inpmat - sqrt(chi^2 + 4 * inpmat * chi * (1 - pi))) / (2 * (chi + s))
  
  LLA <- ifelse(WS.LL < 0, 0, WS.LL)
  ULA <- ifelse(WS.UL > 1, 1, WS.UL)
  
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
