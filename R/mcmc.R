#' MCMC with Metropolis-Hastings
#'
#' A wrapper function to apply the Metropolis-Hastings algorithm for MCMC.
#'
#' @param A A matrix or data frame for your dataset.
#' @param beta_ini Initial parameter values for the MCMC.
#' @param y The vector of observed data.
#' @param iterations Number of MCMC iterations.
#' @param burn_in Number of burn-in iterations to discard.
#' @return An object of class `mcmc` containing the trace of MCMC samples.
#' @export
prog2 <- function(A, beta_ini, y, iterations, burn_in) {
  mcmc_r <- Metro_Hastings(
    li_func = stsu_likelihood,
    pars = beta_ini,
    par_names = c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5"),
    iterations = iterations,
    burn_in = burn_in,
    adapt_par = c(100, 20, 0.5, 0.75),
    y = y,
    nL = nL,
    zA = zA
  )
  return(mcmc(mcmc_r$trace))
}

#' Gelman-Rubin Diagnostics for MCMC Convergence
#'
#' This function performs Gelman-Rubin diagnostics on two MCMC chains.
#' It provides convergence diagnostics using the Gelman-Rubin statistic.
#'
#' @param chain1 MCMC trace from the first chain.
#' @param chain2 MCMC trace from the second chain.
#' @return A list with diagnostics and plots.
#' @export
gelman_diagnostics <- function(chain1, chain2) {
  combined_chains <- mcmc.list(chain1, chain2)
  diag_result <- gelman.diag(combined_chains)
  gelman_plot <- gelman.plot(combined_chains)

  summary_combined <- summary(combined_chains)

  list(
    diagnostics = diag_result,
    plot = gelman_plot,
    summary = summary_combined
  )
}

#' Example Usage for MCMC Convergence Diagnostics
#'
#' This function runs an example of MCMC and Gelman-Rubin diagnostics.
#'
#' @param A Dataset for running the MCMC.
#' @param y Observed data.
#' @param iterations Number of iterations for MCMC.
#' @param burn_in Burn-in period for MCMC.
#' @return Gelman-Rubin diagnostics and related plots.
#' @export
run_diagnostics <- function(A, y, iterations, burn_in) {
  chain1 <- prog2(A = A, beta_ini = c(0.2, 1, 0.6, 2, 0.5, 3), y = y, iterations = iterations, burn_in = burn_in)
  chain2 <- prog2(A = A, beta_ini = c(0.3, 0.9, 0.5, 2, 0.5, 2.9), y = y, iterations = iterations, burn_in = burn_in)

  # Perform Gelman-Rubin diagnostics
  result <- gelman_diagnostics(chain1, chain2)

  # Return result
  result
}
