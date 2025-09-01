#' Likelihood for STSU Weibul Model
#'
#' This function computes the likelihood for the STSU Weibul model using the provided parameters.
#'
#' @param pars A vector of parameters for the Weibul model.
#' @param y A vector of observed event times.
#' @param nL A vector of partition sizes.
#' @param zA A vector of normalized partition weights.
#' @return The likelihood value for the STSU Weibul model.
#' @export
stsu_likelihood <- function(pars, y, nL, zA) {
  f <- function(t1) {
    mu * lambda[k] * t1^(mu - 1) * exp(-t1^mu * lambda[k])
  }
  lambda <- exp(M %*% pars)
  lambda[lambda > exp(700)] <- exp(700) # To avoid overflow
  L <- 10
  LAi <- rep(0, 1)
  LAl <- rep(0, 1)
  int <- rep(0, 1)
  k <- 0
  for (l in 1:L) {
    for (i in 1:nL[l]) {
      k <- k + 1
      int[k] <- integrate(f, 0, y[k], stop.on.error = FALSE)$value
      LAi[i] <- zA[l] * (mu * lambda[k] * (y[k])^(mu - 1) * exp(-(y[k])^mu * lambda[k]))^(nu[k]) *
        (1 - int[k])^(1 - nu[k])
    }
    LAl[l] <- sum(log(LAi))
  }
  LA <- sum(LAl)
  return(LA + prior_beta(pars))
}

#' Log Likelihood for STSU Weibul Model
#'
#' This function computes the log-likelihood for the STSU Weibul model.
#'
#' @param pars A vector of parameters for the Weibul model.
#' @param y A vector of observed event times.
#' @param nL A vector of partition sizes.
#' @param zA A vector of normalized partition weights.
#' @return The log-likelihood value for the STSU Weibul model.
#' @export
loglikelihood <- function(pars, y, nL, zA) {
  f <- function(t1) {
    mu * lambda[k] * t1^(mu - 1) * exp(-t1^mu * lambda[k])
  }
  lambda <- exp(M %*% pars)
  lambda[lambda > exp(700)] <- exp(700) # To avoid overflow
  L <- 10
  LAi <- rep(0, 1)
  LAl <- rep(0, 1)
  int <- rep(0, 1)
  k <- 0
  for (l in 1:L) {
    for (i in 1:nL[l]) {
      k <- k + 1
      int[k] <- integrate(f, 0, y[k], stop.on.error = FALSE)$value
      LAi[i] <- zA[l] * (mu * lambda[k] * (y[k])^(mu - 1) * exp(-(y[k])^mu * lambda[k]))^(nu[k]) *
        (1 - int[k])^(1 - nu[k])
    }
    LAl[l] <- sum(log(LAi))
  }
  LA <- sum(LAl)
  return(LA)
}
