#' BMDU: Multinomial–Dirichlet Unequal Prior Bayesian Method
#'
#' Computes the Bayesian Dirichlet posterior for a multinomial vector
#' using unequal prior parameters. The prior is constructed by dividing
#' the categories into two groups, assigning random priors from different ranges
#' to simulate unequal information across categories.
#'
#' @param x Integer vector of observed counts. Must be non-negative.
#' @param d Integer scalar controlling how the categories are divided into
#'   two groups for constructing unequal Dirichlet priors.
#'
#' @return Prints posterior means, lower and upper 95% credible limits for each category,
#' and the product of the interval widths (volume).
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 2
#' BMDU(y, z)
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats runif quantile
#' @export
BMDU <- function(x, d) {
  k <- length(x)
  
  # Input checks
  if (any(x < 0)) {
    stop("All elements of 'x' must be non-negative integers.")
  }
  if (d < 1 || d > k) {
    stop("The size of the division ('d') should be between 1 and the number of categories.")
  }
  
  # Division for unequal Dirichlet priors
  s1 <- floor(k / d)
  d1 <- runif(s1, 0, 1)       # First subset
  d2 <- runif(k - s1, 1, 2)   # Second subset
  a  <- c(d1, d2)
  
  # Posterior simulation
  p  <- x + a
  dr <- MCMCpack::rdirichlet(10000, p)
  
  # Initialize vectors
  m    <- numeric(k)
  l    <- numeric(k)
  u    <- numeric(k)
  diff <- numeric(k)
  
  for (j in 1:k) {
    l[j]    <- round(quantile(dr[, j], 0.025), 4)
    u[j]    <- round(quantile(dr[, j], 0.975), 4)
    m[j]    <- round(mean(dr[, j]), 4)
    diff[j] <- u[j] - l[j]
  }
  
  # Compute volume
  vol <- prod(diff)
  
  # Combine results into a single data frame
  result <- data.frame(
    mean      = m,
    lower_95  = l,
    upper_95  = u,
    interval  = diff,
    volume    = vol
  )
  
  return(result)
}
