#' Prior Distribution for Beta Parameters
#'
#' This function computes the prior distribution for the beta parameters in a survival model.
#'
#' @param pars A vector of beta parameters
#' @return The sum of the log of the normal prior density for each beta
#' @export
prior_beta <- function(pars) {
  p <- length(pars)
  beta <- pars[1:p]
  prior <- dnorm(beta, 0, 10)
  return(sum(prior))
}
