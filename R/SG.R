#' SG: Sison and Glaz Method for Simultaneous Confidence Intervals
#'
#' Computes simultaneous confidence intervals for multinomial proportions
#' using the Sison and Glaz (1995) method. The function implements the
#' truncated Poisson approach, approximates the required probabilities via
#' Edgeworth expansion, and determines the limits that ensure the overall
#' confidence level (1 - alpha).
#'
#' @param x Integer vector of observed cell counts corresponding to a
#' categorical dataset. All entries must be non-negative.
#' @param alpha Desired statistical significance level.
#'
#' @return Prints the original and adjusted confidence intervals for each
#' category, along with the volume (product of interval widths). 
#'
#' @details
#' This function implements the simultaneous confidence interval construction
#' proposed by Sison and Glaz (1995). It is based on a truncated Poisson model
#' with factorial and central moment calculations, Edgeworth expansion for
#' probability approximation, and adjustment of limits to ensure they remain
#' within the [0,1] range.
#'
#' The computed volume represents the product of the widths of all confidence
#' intervals and serves as a measure of the overall uncertainty.
#'
#' @references
#' Sison, C. P., and Glaz, J. (1995).
#' *Simultaneous Confidence Intervals and Sample Size Determination for Multinomial Proportions.*
#' Journal of the American Statistical Association, **90**, 366–369.
#'
#' @author
#' Dr. M. Subbiah
#'
#' @seealso
#' \code{\link{BMDE}}, \code{\link{WALD}}, \code{\link{GM}}
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 0.05
#' SG(y, z)
#'
#' @importFrom stats ppois dpois
#' @export
SG <- function(x, alpha) {
  if (any(x < 0)) stop("All elements of 'x' must be non-negative.")
  if (sum(x) == 0) stop("Sum of 'x' must be positive.")
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be between 0 and 1.")
  
  sgp <- function(c) {
    s <- sum(x)
    k <- length(x)
    b <- x - c
    a <- x + c
    
    ## Factorial moments for truncated Poisson
    fm1 <- fm2 <- fm3 <- fm4 <- numeric(k)
    for (i in 1:k) {
      denom <- (ppois(a[i], x[i]) - ppois(b[i] - 1, x[i]))
      fm1[i] <- x[i] * (ppois(a[i] - 1, x[i]) - ppois(b[i] - 2, x[i])) / denom
      fm2[i] <- x[i]^2 * (ppois(a[i] - 2, x[i]) - ppois(b[i] - 3, x[i])) / denom
      fm3[i] <- x[i]^3 * (ppois(a[i] - 3, x[i]) - ppois(b[i] - 4, x[i])) / denom
      fm4[i] <- x[i]^4 * (ppois(a[i] - 4, x[i]) - ppois(b[i] - 5, x[i])) / denom
    }
    
    ## Central moments for truncated Poisson
    m1 <- fm1
    m2 <- fm2 + fm1 - fm1^2
    m3 <- fm3 + fm2 * (3 - 3 * fm1) + (fm1 - 3 * fm1^2 + 2 * fm1^3)
    m4 <- fm4 + fm3 * (6 - 4 * fm1) +
      fm2 * (7 - 12 * fm1 + 6 * fm1^2) +
      fm1 - 4 * fm1^2 + 6 * fm1^3 - 3 * fm1^4
    m4t <- m4 - 3 * m2^2
    
    s1 <- sum(m1)
    s2 <- sum(m2)
    s3 <- sum(m3)
    s4 <- sum(m4t)
    
    ## Edgeworth expansion
    g1 <- s3 / (s2^(3 / 2))
    g2 <- s4 / (s2^2)
    z <- (s - s1) / sqrt(s2)
    z2 <- z^2; z3 <- z^3; z4 <- z^4; z6 <- z^6
    poly <- 1 + g1 * (z3 - 3 * z) / 6 + g2 * (z4 - 6 * z2 + 3) / 24 +
      (g1^2) * (z6 - 15 * z4 + 45 * z2 - 15) / 72
    f <- poly * exp(-z2 / 2) / sqrt(2 * pi)
    
    ## Probability computation
    pc <- numeric(k)
    for (i in 1:k) {
      pc[i] <- ppois(a[i], x[i]) - ppois(b[i] - 1, x[i])
    }
    
    pcp <- prod(pc)
    pps <- 1 / dpois(s, s)
    rp <- pps * pcp * f / sqrt(s2)
    rp
  }
  
  s <- sum(x)
  c <- 1:s
  y <- sapply(c, function(ci) round(sgp(ci), 4))
  
  vc <- which(diff(y < (1 - alpha)) != 0)[1]
  if (is.na(vc)) stop("Could not determine critical value 'c'. Try different alpha or data.")
  delta <- ((1 - alpha) - y[vc]) / (y[vc + 1] - y[vc])
  
  sp <- x / s
  LL <- round(sp - (vc / s), 4)
  UL <- round(sp + (vc / s) + (2 * delta / s), 4)
  
  LLA <- ifelse(LL < 0, 0, LL)
  ULA <- ifelse(UL > 1, 1, UL)
  
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
