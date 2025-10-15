#' @importFrom MCMCpack rdirichlet
NULL

#' BMDE: Multinomial–Dirichlet Equal Prior Bayesian Method
#'
#' Computes the Bayesian Dirichlet posterior for a multinomial vector
#' with equal prior parameters and returns the posterior mean, 95% credible intervals,
#' and the volume of those intervals.
#'
#' @param x Integer vector of observed counts. Must be non-negative.
#' @param p Numeric scalar or vector specifying Dirichlet prior parameters. Must be non-negative.
#'
#' @return Prints posterior means, lower and upper 95% credible limits for each category,
#' and the product of the interval widths (volume).
#'
#' @examples
#' y <- c(44, 55, 43, 32, 67, 78)
#' z <- 1
#' BMDE(y, z)
#'
#' @importFrom MCMCpack rdirichlet
#' @export
BMDE <- function(x, p) {
  k <- length(x)
  n_r <- 10000
  
  # Input check
  if (any(x < 0)) {
    stop("All elements of 'x' must be non-negative integers.")
  }
  
  # Dirichlet posterior simulation
  po <- x + p
  dr <- MCMCpack::rdirichlet(n_r, po)
  
  # Initialize result vectors
  mean_post <- numeric(k)
  lower     <- numeric(k)
  upper     <- numeric(k)
  interval  <- numeric(k)
  
  # Compute posterior summaries
  for (j in 1:k) {
    mean_post[j] <- round(mean(dr[, j]), 4)
    lower[j]     <- round(quantile(dr[, j], 0.025), 4)
    upper[j]     <- round(quantile(dr[, j], 0.975), 4)
    interval[j]  <- upper[j] - lower[j]
  }
  
  # Combine results into a single data frame
  res <- data.frame(
    mean      = mean_post,
    lower_95  = lower,
    upper_95  = upper,
    interval  = interval,
    volume = prod(interval)
  )
  
  return(res)
}
