#' Get re-fitted scores based on FPCA basis.
#'
#' Get scores based on a linear mixed effects model using the FPCA
#' eigenfunctions as a basis for time.
#'
#' @param fpca Object from FPCA fitting.
#' @param df Number of degrees of freedom to use.
#'
#' @return A data frame containing the scores.
#'
#' @export
#'

fpme <- function(fpca, ids = NULL, df) {
  yds <- get_ymu(fpca, ids = ids)
  phi <- get_phi_i(fpca, yds$tt)
  fit <- lme4::lmer(ymu ~ -1 + phi[,1:df] + (-1 + phi[,1:df] | id), data = yds)
  b <- lme4::ranef(fit)$id
  colnames(b) <- NULL
  data.frame(id = rownames(b),
             xime = b)
}
