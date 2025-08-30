#' Summary method for GCElm objects
#'
#' @description
#' Provides a summary of a model fitted with \code{GCElm}, including estimated
#' coefficients, standard errors, test statistics, and p-values.
#'
#' @param object An object of class \code{"GCElm"} returned by \code{GCElm}.
#' @param correlation logical; if TRUE, the correlation matrix of the estimated parameters is returned and printed.
#' @param symbolic.cor logical; If TRUE, print the correlations in a symbolic form (see symnum) rather than as numbers.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' A \code{"summary.GCElm"} object with components described below:
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{coefficients}}{A matrix with estimates, standard errors, t values, and p-values.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1234)
#' N <- 1000; K <- 10
#' betas <- c(-1, 1, -1, 1, rep(0, K-4))
#' X <- matrix(runif(N*K), nrow = N, ncol = K)
#' y <- as.vector(X %*% betas + rnorm(N))
#' Z <- matrix(rep(seq(-10, 10, 2.5), K + 1), nrow = K + 1, byrow = TRUE)
#' df <- data.frame(y, X)
#' fit <- GCElm(y ~ ., data = df, Z = Z)
#' }
#'
#' @importFrom stats var
#' @importFrom stats pt
#' @importFrom stats coef
#' @method summary GCElm
#' @export
summary.GCElm <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.GCElm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL, 
                                                                 c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
    if (correlation) 
      ans$correlation <- ans$cov.unscaled
    return(ans)
  }
  if (is.null(z$terms)) 
    stop("invalid 'GCElm' object:  no 'terms' component")
  if (!inherits(object, "GCElm")) 
    warning("calling summary.GCElm(<fake-GCElm-object>) ...")
  Qr <- object$qr
  n <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not a \"GCElm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  if (!is.null(z$offset)) {
    f <- f - z$offset
  }
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) 
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) * 
      1e-30) 
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  #se <- sqrt(diag(R) * resvar)
  se <- sqrt(diag(z$var_beta))/z$nu
  est <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <- cbind(Estimate = est, `Std. Error` = se, 
                            `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf, 
                                                                  lower.tail = FALSE))
  ans$aliased <- is.na(z$coefficients)
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, NCOL(Qr$qr))
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
                                                       df.int)/rdf)
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
                        numdf = p - df.int, dendf = rdf)
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
                                                             1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action)) 
    ans$na.action <- z$na.action
  class(ans) <- c("summary.GCElm", "summary.lm")
  ans
}

