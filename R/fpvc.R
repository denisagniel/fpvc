#' Perform testing of summarized longitudinal trajectory.
#'
#' Perform testing of summarized longitudinal trajectory.
#'
#' @param data Data frame containing scores summarizing longitudinal trajectories,
#'     genetic markers, and any covariates to adjust for, one row per
#'     individual.
#' @param xicols Character vector indicating columns of \code{data} that
#'     correspond to the scores summarizing longitudinal trajectory.
#' @param zcols Character vector indicating columns of \code{data} that
#'     correspond to the genetic markers.
#' @param xcols Character vector indicating columns of \code{data} that
#'     correspond to the covariates to adjust for.
#'
#' @return A p-value.
#'
#' @export
#'
#'
fpvc <- function(data, xicols, zcols, xcols = NULL, model) {
  xi <- dplyr::select_(data, xicols)
  xi <- as.matrix(xi)
  z <- dplyr::select_(data, zcols)
  z <- as.matrix(z)
  if (!is.null(xcols)) {
    x <- dplyr::select_(data, xcols)
    x <- as.matrix(x)
  } else x <- NULL
  scoreTest(y = xi, z = z, x = x, model = model)
}
