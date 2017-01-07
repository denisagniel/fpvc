#' Generate data for simulation.
#'
#' Generate data with varying amounts of nonlinearity and complexity consistent
#' with the simulations for "Functional principle variance component testing
#' for a genetic association study of HIV progression."
#'
#' @param n Number of samples.
#' @param nt Average number of time points.
#' @param maf Minor allele frequency of SNP.
#' @param beta Effect size.
#' @param TT Length of study period.
#' @param alpha Parameter controlling the linearity of genetic effect, between 0
#' (nonlinear effect) and 1 (linear effect).
#' @param beta Parameter controlling the complexity of the mean process between
#' 0 (simple mean process) and 1 (complex mean process)
#'
#' @importFrom stats rnorm
#'
#' @export
#'

generateData <- function(n, nt, maf, beta, TT, alpha, gamma) {
  g = sample(0:2, prob = c((1-maf)^2, 2*maf*(1-maf), maf^2), size = n, replace = TRUE)
  b = rnorm(n)
  b2 = rnorm(n)
  b3 = rnorm(n)
  n_t = stats::rpois(n, nt-2) + 2

  sigma = 0.5
  sigma1 = 0.5
  sigma2 = 0.5
  y1 <- y2 <- y1ds <- y2ds <- tt <- list()
  for (i in 1:n) {
    thisn <- n_t[i]
    thisg <- g[i]
    tt[[i]] <- sort(stats::runif(thisn)*TT)
    y1[[i]] <- sin(tt[[i]]) + gamma*(sin(tt[[i]]/3) + cos(tt[[i]])) +
      (1-gamma)*(b2[i]*cos(tt[[i]]/4)*sigma1 + b[i]) +
      g[i]*beta*(alpha*(cos(tt[[i]]) + cos(tt[[i]]/10) - sin(3*tt[[i]])) +
                   (1-alpha)*tt[[i]]/7) + rnorm(thisn)*sigma
    y1ds[[i]] <- data.frame(y1 = y1[[i]],
      g = g[i],
      id = i,
      tt = tt[[i]])
    y2[[i]] = sin(tt[[i]]) + gamma*(-sin(tt[[i]]/3) - cos(tt[[i]])) +
      (1-gamma)*(b3[i]*cos(tt[[i]]/4)*sigma2 + b[i]) -
      g[i]*beta*(alpha*(cos(tt[[i]]) + cos(tt[[i]]/10) - sin(3*tt[[i]])) +
                   (1-alpha)*tt[[i]]/7) + rnorm(thisn)*sigma
      y2ds[[i]] <- data.frame(y2 = y2[[i]],
                         g = g[i],
      id = i,
      tt = tt[[i]])
  }
  mu1 <- mean(unlist(y1))
  sd1 <- stats::sd(unlist(y1))
  mu2 <- mean(unlist(y2))
  sd2 <- stats::sd(unlist(y2))
  y1_scaled <- lapply(y1, function(x) (x - mu1)/sd1)
  y2_scaled <- lapply(y2, function(x) (x - mu2)/sd2)
  ds1 <- dplyr::bind_rows(y1ds)
  ds2 <- dplyr::bind_rows(y2ds)
  # list(y1 = y1, y2 = y2, y1_scaled = y1_scaled,
  #      y2_scaled = y2_scaled,
  #      tt = tt, ds = merge(ds1, ds2))
  merge(ds1, ds2)
}

