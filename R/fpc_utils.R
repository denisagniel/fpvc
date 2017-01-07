get_mu_i <- function(fpca, tt = NULL) {
  if (is.null(tt)) {
    tt <- fpca$inputData$Lt %>% unlist
  }
  stats::approx(fpca$workGrid, fpca$mu, xout = tt, rule = 2)$y
}



listdata_to_ds <- function(fpca, ids = NULL) {
  yy <- unlist(fpca$inputData$Ly)
  tt <- unlist(fpca$inputData$Lt)
  n <- length(fpca$inputData$Ly)
  n_i <- unlist(lapply(fpca$inputData$Ly, length))
  if (is.null(ids)) {
    id <- rep(1:n, n_i)
  } else id <- rep(ids, n_i)
  data.frame(id = id, tt = tt, yy =yy)
}



get_ymu <- function(fpca, ids = NULL, y = NULL, tt = NULL) {
  if (is.null(y) | is.null(tt)) {
    yds <- listdata_to_ds(fpca, ids = ids)
    dplyr::mutate(yds, mu = get_mu_i(fpca, tt),
                   ymu = yy - mu)
  } else {
    mui <- get_mu_i(fpca, tt)
    y - mui
  }
}

