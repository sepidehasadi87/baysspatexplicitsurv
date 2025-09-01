#' Metropolis-Hastings MCMC for STSU Weibul Model
#'
#' This function runs the Metropolis-Hastings MCMC algorithm for fitting the STSU Weibul model.
#'
#' @param A A matrix or data frame for your dataset.
#' @param beta_ini Initial parameter values for the MCMC.
#' @param y The vector of observed data.
#' @param iterations Number of MCMC iterations.
#' @param burn_in Number of burn-in iterations to discard.
#' @return An object of class `mcmc` containing the trace of MCMC samples.
#' @export
prog2 <- function(A, beta_ini, y, iterations, burn_in) {
  # Running MCMC using the Metro-Hastings algorithm
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

  # Calculate parameter estimates and standard deviations
  beta_hat <- apply(mcmc_r$trace, 2, mean)
  SD <- apply(mcmc_r$trace, 2, sd)
  print(round(beta_hat, 3))
  print(round(SD, 3))

  return(mcmc_r)
}

#' Compute the Mean Squared Error (MSE)
#'
#' This function computes the Mean Squared Error (MSE) between estimated and true parameter values.
#'
#' @param estimate Estimated parameter value.
#' @param true_value True parameter value.
#' @return The computed MSE.
#' @export
MSE <- function(estimate, true_value) {
  return((estimate - true_value)^2)
}

#' Run MCMC Loop for Model Fitting
#'
#' This function runs the MCMC loop for a specified number of iterations and calculates various diagnostics (DIC, BIC, etc.).
#'
#' @param A The data matrix for the model.
#' @param y The observed data.
#' @param iterations The number of MCMC iterations.
#' @param burn_in The number of burn-in iterations.
#' @return A list of calculated metrics (DIC, BIC, pD, etc.).
#' @export
run_mcmc_loop <- function(A, y, iterations = 10000, burn_in = 1000, a = 50) {
  # Initialize storage for results
  dicS <- rep(0, a)
  betaS <- matrix(0, a, 6)
  bicS <- rep(0, a)

  D_bar <- rep(0, a)
  pD2 <- rep(0, a)
  DIC2 <- rep(0, a)

  # True values of parameters
  beta_real <- c(1, 1.2, 0.5, 3, 0.2, 5)

  # Main MCMC loop
  for (i in 1:a) {
    cat("Iteration No. =", i, "\n")

    # Run MCMC for the current iteration
    mcmc_r <- prog2(A = A, beta_ini = c(0.2, 1, 0.6, 2, 0.5, 3), y = y, iterations = iterations, burn_in = burn_in)
    mcmc_r1 <- mcmc_thin(mcmc_r)

    # Store BIC and DIC values
    bicS[i] <- BCI(mcmc_r1)
    betaS[i, ] <- apply(mcmc_r1$trace, 2, mean)
    dicS[i] <- mcmc_r1$DIC

    # Compute pD and DIC
    D_bar[i] <- (-2) * loglikelihood(pars = betaS[i, ], y = y, nL = nL, zA = zA)
    var <- 0
    G <- 8835
    for (g in 1:G) {
      var <- var + ((-2) * loglikelihood(pars = mcmc_r$trace[g, ], y = y, nL = nL, zA = zA) - D_bar[i])^2
    }
    pD2[i] <- (1 / 2) * (1 / (G - 1)) * var
    DIC2[i] <- D_bar[i] + 2 * pD2[i]
  }

  # Calculate summary metrics
  pDhat <- round(mean(pD2), 3)
  DIC2 <- round(mean(DIC2), 3)
  DICS <- round(mean(dicS), 3)

  # Calculate BIC and Betahat
  BICt <- sum(bicS)
  Betahat <- round(BICt[, 3] / a, 3)
  BIC <- round(BICt / a, 3)

  # Calculate MSE for all parameters
  MSES <- sapply(1:6, function(i) MSE(Betahat[i], beta_real[i]))
  MSES <- round(MSES, 3)

  # Return results
  return(list(DIC = DICS, Betahat = Betahat, BIC = BIC, MSE = MSES))
}

#' Thin the MCMC Trace
#'
#' This function thins the MCMC trace by selecting one sample for every specified interval.
#'
#' @param mcmc_r The MCMC result object.
#' @return A thinned MCMC trace.
#' @export
mcmc_thin <- function(mcmc_r) {
  # Thinning the MCMC trace (keep every nth iteration)
  thin_trace <- mcmc_r$trace[seq(1, nrow(mcmc_r$trace), by = 10), ]
  return(list(trace = thin_trace, DIC = mcmc_r$DIC))
}
