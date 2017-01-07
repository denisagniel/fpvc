#' Variance component score test.
#'
#' Tests for an association between multivariate \code{y} and multivariate \code{z} conditional
#' on multivariate \code{x}.
#'
#' @param y The outcomes of interest, n x m matrix.
#' @param z The predictors of interest, n x p matrix.
#' @param x Covariates to adjust for, n x q matrix.
#' @param model Model used to adjust for \code{x}, currently only proportional
#' odds ('po'), binomial, logistic, and probit supported.
#'
#' @importFrom stats predict glm
#'
#' @return p-value of the test
#'
#' @examples
#' n <- 1000
#' p <- 2
#' q <- 5
#' m <- 3
#' x <- rnorm(n*q) %>% matrix(n, q)
#' pi <- pnorm(x)
#' z <- rbinom(n*p, size = 1, prob = pi) %>% matrix(n, p)
#' mu <- rowMeans(x) + z[,1]/10
#' y <- rnorm(n*m, mean = mu) %>% matrix(n, m)
#' scoreTest(y = y, z = z, x = x, model = 'probit')
#'
#' @export
#'
scoreTest <- function(y, z, x = NULL, model) {
  y <- as.matrix(y)
  z <- as.matrix(z)
  n <- nrow(z)
  if (!is.null(x)) {
    x <- as.matrix(x)
  } else {
    x <- matrix(0, n, 0)
  }
  x1 <- cbind(1, x)
  p <- ncol(z)
  q <- ncol(x) + 2
  m <- ncol(y)
  stopifnot(nrow(y) == n)
  stopifnot(nrow(x) == n)

  ## dgb is n x q, this is the first derivative of link function for E[z | x]
  ## A.xx is the information matrix
  ## U is the score function, s.t. A.xx^-1 x U is the influence fn
  S <- lapply(1:p, function(j) {
    if (model == 'po') {
      fit.x <- rms::lrm(z[,j] ~ x1 - 1)
      dgb <- dg(x, fit.x, model = 'po') ## n x q
      A.xx <- fit.x$info.matrix
    } else if (model == 'binomial') {
      fit.x <- stats::glm(cbind(z[,j], 2-z[,j]) ~ x1 -1, family = binomial)
      dgb <- dg(x, fit.x, model = 'binomial') ## n x q
      A.xx <- t(x1)%*%(dgb)
    } else if (model == 'logistic') {
      fit.x <- stats::glm(z[,j] ~ x1 -1, family = binomial)
      dgb <- dg(x, fit.x, model = 'logistic')
      # A.xx <- fit.x$info.matrix
      A.xx <- t(x1)%*%(dgb)
    } else if (model == 'probit') {
      # browser()
      fit.x <- stats::glm(z[,j] ~ x1 -1, family = stats::binomial(link = 'probit'))
      dgb <- dg(x, fit.x, model = 'probit')
      A.xx <- t(x1) %*% diag(fit.x$weights) %*% x1
    } else stop("Please choose one of {'po', 'binomial', 'logistic', 'probit'}")
    res.x <- z[,j] - g(fit.x, model = model) ## n x 1
    # browser()
    R.yx <- t(y) %*% (dgb) ## m x q
    S.hat <- colSums(y*res.x)/sqrt(n)
    S.exp <- sqrt(n/(n-q)) * ((y*res.x) - t(R.yx %*% solve(A.xx) %*%
                                              t(U(x, z[,j], fit.x, model = model))))
    list(S.hat = S.hat, S.exp = S.exp)
  })
  S.hat <- unlist(lapply(S, function(x) x$S.hat))
  S.exp <- do.call(cbind, lapply(S, function(x) x$S.exp))
  Sig = var(S.exp); qhat = sum(S.hat^2);
  tmpsvd = eigen(Sig,symmetric=T);
  lam.all = tmpsvd$val
  CompQuadForm::davies(qhat,lam.all)$Qq
}

expit <- function(xx){exp(xx)/(1+exp(xx))}
d.expit <- function(xx){expit(xx)*(1-expit(xx))}

g <- function(fit, model) {
  if (model == 'po') {
    fitted.p <- predict(fit, type = 'fitted')
    rowSums(fitted.p)
  } else if (model == 'binomial') {
    2*predict(fit, type = 'response')
  } else if (model %in% c('logistic', 'probit')) {
    predict(fit, type = 'response')
  } else stop("Please choose one of {'po', 'binomial', 'logistic', 'probit'}")
}


dg <- function(xx, fit, model) {
  if (model == 'po') {
    a1xb <- predict(fit)
    a2xb <- predict(fit, kint = 2)
    matrix(c(d.expit(a1xb), d.expit(a2xb),
             (d.expit(a1xb) + d.expit(a2xb)) * xx ), nrow(xx), 2+ncol(xx))
  } else if (model %in% c('binomial', 'logistic')) {
    xb <- predict(fit)
    ifelse(model == 'binomial', 2, 1)*d.expit(xb) * cbind(1, xx) ## if binomial, then g = 2*pi, so the derivative has a factor of 2
  } else if (model == 'probit') {
    xb <- predict(fit)
    dnorm(xb)*cbind(1, xx)
  } else stop("Please choose one of {'po', 'binomial', 'logistic', 'probit'}")
}

U <- function(xx, zz, fit, model) {
  # browser()
  if (model == 'po') {
    p.fit <- predict(fit, type = 'fitted')
    p.0 <- 1 - p.fit[,1]
    p.1 <- p.fit[,1] - p.fit[,2]
    p.2 <- p.fit[,2]
    a1xb <- predict(fit)
    a2xb <- predict(fit, kint = 2)
    z.0 <- 1*(zz == 0)
    z.1 <- 1*(zz == 1)
    z.2 <- 1*(zz == 2)
    U.a1 <- -d.expit(a1xb)*(z.0/p.0 - z.1/p.1)
    U.a2 <- -d.expit(a2xb)*(z.1/p.1 - z.2/p.2)
    U.b <- z.0/p.0*(-d.expit(a1xb))*xx +
      z.1/p.1*(d.expit(a1xb) - d.expit(a2xb))*xx +
      z.2/p.2*(d.expit(a2xb))*xx
    cbind(U.a1, U.a2, U.b)
  } else if (model == 'binomial') {
    # browser()
    matrix(zz - 2*predict(fit, type = 'response'), nrow(xx), ncol(xx) + 1) * cbind(1, xx)
  } else if (model == 'logistic') {
    matrix(zz - predict(fit, type = 'response'),
           nrow(xx), ncol(xx) + 1) * cbind(1, xx)
  } else if (model == 'probit') {
    pi_i <- predict(fit, type = 'response')
    xb <- predict(fit)
    matrix(stats::dnorm(xb)*(zz - pi_i)/pi_i/(1-pi_i),
           nrow(xx), ncol(xx) + 1) * cbind(1, xx)
  } else stop("Please choose one of {'po', 'binomial', 'logistic', 'probit'}")
}
