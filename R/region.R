#' Generate a Region
#'
#' This function generates a region of random points and plots them.
#'
#' @param n Number of points to generate.
#' @param seed Seed for random number generation.
#' @return A matrix of points.
#' @export
region <- function(n, seed) {
  set.seed(seed)
  A <- cbind(c(runif(n, min = 0, max = 1)), c(runif(n, min = 0, max = 1)))
  plot(A)
  abline(v = 0.2, col = "red")
  abline(v = 0.4, col = "red")
  abline(v = 0.6, col = "red")
  abline(v = 0.8, col = "red")
  abline(h = 0.5, col = "red")
  A
}

# R/data.R

#' Simulate Data for Survival Analysis
#'
#' This function generates simulated data for survival analysis, including covariates,
#' time-to-event, and censoring times.
#'
#' @param n Number of observations
#' @param A A matrix representing the region for spatial covariates
#' @return A matrix of simulated data with covariates, time-to-event, and censoring information
#' @export
data <- function(n, A) {
  library(geoR)

  L <- 10 # number of partitions
  p <- 6 # number of covariates + 1
  ng <- 9 # number of grids in each partition

  # Generate p covariates m1,..., mp
  m1 <- rbinom(n, 1, 0.5)
  m2 <- rbinom(n, 1, 0.5)
  m3 <- rnorm(n, 0, 1)
  m5 <- grf(n, grid = A, cov.model = "gaussian", cov.pars = c(sigmasq = 4, phi = 0.25, tausq = 0))
  m4 <- grf(n, grid = A, cov.model = "gaussian", cov.pars = c(sigmasq = 10, phi = 0.25, tausq = 0))

  one <- rep(1, n)
  M <- matrix(c(one, m1, m2, m3, m4$data, m5$data), nrow = n, byrow = FALSE)

  # Generate time to event
  beta_real <- c(1, 1.2, 0.5, 3, 0.2, 5)
  lambda <- exp(M %*% beta_real)
  tt <- rexp(n = n, rate = lambda)

  # Generate censoring time
  cc <- runif(n, 0, 2)

  # Generate y and nu
  y <- rep(0, n)
  nu <- rep(0, n)

  for (i in 1:n) {
    y[i] <- min(cc[i], tt[i])
    if (tt[i] <= cc[i]) {
      nu[i] <- 1
    } else {
      nu[i] <- 0
    }
  }

  d <- cbind(M, y, nu, tt, cc)
  return(d)
}
