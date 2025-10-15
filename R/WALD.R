#' WALD: Wald Method for Simultaneous Confidence Intervals
#'
#' Computes simple Wald-type simultaneous confidence intervals for multinomial
#' proportions. These intervals are symmetric about the sample proportions and
#' do not use continuity corrections, thus avoiding zero-width intervals even
#' for extreme sample proportions.
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
#' The adjusted limits are truncated to stay within the [0,1] range.
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
#' \code{\link{BMDE}}, \code{\link{WALDCC}}, \code{\link{SG}}
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 0.05
#' WALD(y, z)
#'
#' @importFrom stats qchisq
#' @export
WALD <- function(inpmat, alpha) {
  if (any(inpmat < 0)) stop("All elements of 'inpmat' must be non-negative.")
  if (sum(inpmat) == 0) stop("Sum of 'inpmat' must be positive.")
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be between 0 and 1.")
  
  k <- length(inpmat)
  s <- sum(inpmat)
  chi <- qchisq(1 - alpha, df = 1)
  pi <- inpmat / s
  
  WALD.LL <- pi - sqrt(chi * pi * (1 - pi) / s)
  WALD.UL <- pi + sqrt(chi * pi * (1 - pi) / s)
  
  LLA <- ifelse(WALD.LL < 0, 0, WALD.LL)
  ULA <- ifelse(WALD.UL > 1, 1, WALD.UL)
  
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
