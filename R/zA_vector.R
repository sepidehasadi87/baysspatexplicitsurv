#' Calculate zA and nL Vectors
#'
#' This function calculates the `zA` vector and partition sizes (`nL`) based on input data `A`.
#'
#' @param n The number of data points.
#' @param A A matrix or data frame representing spatial data points.
#' @return A list containing the partition sizes `nL` and the normalized `zA` vector.
#' @export
zA_vector <- function(n, A) {
  L <- 10
  ng <- 9
  ww <- partition(A) # Assuming partition is another function in your package
  nL <- ww[[1]]
  w <- ww[[2]]
  x <- rnorm(ng, 0, 25)
  zA <- rep(0, L)
  nL1 <- c(0)
  zAT <- 0

  for (l in 1:L) {
    wl <- w[((l - 1) * 9 + 1):(l * 9), ]
    nL1[l + 1] <- nL1[l] + nL[l]
    Al <- A[(nL1[l] + 1):nL1[l + 1], ]
    for (j in 1:ng) {
      dd <- 0
      for (i in 1:9) {
        for (m in 1:length(Al[, 1])) {
          dd <- dd + dist(cbind(c(wl[i, 1], Al[m, 1]), c(wl[i, 2], Al[m, 2])))
        }
      }
      zA[l] <- zA[l] + exp(x[j]) * exp(-dd / 2)
    }
  }
  return(list(nL, zA / sum(zA)))
}
