#' Fit functional principal components analysis
#'
#' Fit FPCA and return scores matched to IDs and return eigenfunctions
#' matched to IDs.
#'
#' @param ds Data frame containing data for analysis, one row per observation.
#' @param ycols Character vector identifying columns of \code{ds} corresponding to longitudinal outcomes.
#' @param tcol Character vector identifying the column of \code{ds} corresponding to time.
#' @param idcol Character vector identifying the column of \code{ds} corresponding to patient ID.
#' @param options List of options to pass to \code{FPCA}, see documentation at \code{\link[fdapace]{FPCA}}.
#'
#' @importFrom dplyr select_ bind_cols
#'
#' @return A list containing a data frame of FPCA scores and IDs, a data frame
#'    containing eigenfunction evaluations at observed time points and IDs,
#'    and a list of FPCA fits.
#'
#' @export
#'

fpca <- function(ds, ycols, tcol, idcol, options) {
  # browser()
  list_data <- dataframe_to_list(ds, ycols, tcol, idcol)
  fpca_l <- with(list_data, runFPCA_list(y_lists, t_list, num_y, options))
  score_l <- list()
  for (ii in 1:length(fpca_l)) {
    score_l[[ii]] <- data.frame(fpca_l[[ii]]$xiEst)
    colnames(score_l[[ii]]) <- stringr::str_c('xi', ii, '.', 1:ncol(score_l[[ii]]))
  }
    # lapply(fpca_l, function(x) data.frame(xi = x$xiEst))
  score_ds <- bind_cols(score_l)
  score_ds$id <- list_data$id

  tt <- unlist(select_(ds, tcol))
  phi_l <- list()
  for (ii in 1:length(fpca_l)) {
    phi_l[[ii]] <- data.frame(get_phi_i(fpca_l[[ii]], tt))
    colnames(phi_l[[ii]]) <- stringr::str_c('phi', ii, '.', 1:ncol(phi_l[[ii]]))
  }
  # phi_l <- lapply(fpca_l, function(fpc) data.frame(phi = get_phi_i(fpc, tt)))
  phi_ds <- bind_cols(phi_l)
  phi_ds$id <- select_(ds, idcol)
  phi_ds$tt <- select_(ds, tcol)

  list(scores = score_ds, phi_i = phi_ds, ids = list_data$id, fpca_fits = fpca_l)
}

dataframe_to_list <- function(ds, ycols, tcol, idcol) {
  y_lists <- lapply(ycols, function(y) {
    plyr::dlply(ds, .variables = idcol, function(x) {
      unlist(dplyr::select_(x, y))
    })
  })
  t_list <- plyr::dlply(ds, .variables = idcol, function(x) {
    unlist(dplyr::select_(x, tcol))
  })
  list(y_lists = y_lists, t_list = t_list, id = names(t_list), num_y = length(ycols))
}

runFPCA_list <- function(y_lists, t_list, num_y, options) {
  if (num_y == 1) {
    fdapace::FPCA(y_lists, t_list, options)
  } else {
    # browser()
    fpca_l <- lapply(y_lists, function(yy) {
      fdapace::FPCA(yy, t_list, options)
    })
    fpca_l
  }
}

get_phi_i <- function(fpca, tt = NULL) {
  if (is.null(tt)) {
    tt <- unlist(fpca$inputData$Lt)
  }
  apply(fpca$phi, 2, function(x) {
    stats::approx(fpca$workGrid, x, xout = tt, rule = 2)$y
  })
}

