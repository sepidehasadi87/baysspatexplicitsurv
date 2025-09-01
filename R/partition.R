#' Partition Spatial Data
#'
#' This function divides a given spatial region into smaller partitions.
#'
#' @param A A matrix representing the region for spatial covariates
#' @return A list containing partition sizes and grid locations
#' @export
partition <- function(A) {
  L <- 10
  nL <- rep(0, 1)
  y1 <- 0
  y2 <- 0
  for (i in 1:2) {
    y1 <- y2
    y2 <- y1 + 0.5
    x1 <- 0
    x2 <- 0
    for (j in 1:5) {
      Al <- c()
      x1 <- x2
      x2 <- x1 + 0.2
      index <- which((A[, 1] > x1) & ((A[, 1] == x2 | A[, 1] < x2)) & (A[, 2] > y1) & ((A[, 2] == y2 | A[, 2] < y2)))
      m <- j + 5 * (i - 1)
      Al <- A[index, ]
      nL[m] <- dim(Al)[1]
    }
  }

  # Generating Partitions Al, nL
  nL1 <- c()
  nL1[1] <- 0
  for (l in 1:L) {
    nL1[l + 1] <- nL1[l] + nL[l]
    Al <- A[(nL1[l] + 1):nL1[l + 1], ]
  }

  # Generating grids w
  w <- c()
  s1 <- matrix(c(0, .1, .2, 0, .1, .2, 0, .1, .2, 0, 0, 0, .25, .25, .25, .5, .5, .5), 9, 2)
  for (i in 1:2) {
    for (j in 1:5) {
      s <- cbind(s1[, 1] + .2 * (j - 1), s1[, 2] + .5 * (i - 1))
      w <- rbind(w, s)
    }
  }
  return(list(nL, w))
}
