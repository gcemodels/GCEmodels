#' Summary of a GCE/GME Binary Logit Model
#'
#' Computes Wald-type inference for the regression coefficients of a
#' GCE/GME binary logit model, including standard errors, z-statistics,
#' and p-values. Inference is based on either a sandwich or a classical
#' covariance estimator.
#'
#' @param object An object of class \code{"GCElogit"}, as returned by
#'   \code{GCElogit()}.
#' @param method Character string specifying the covariance estimator to use.
#'   Must be one of \code{"sandwich"} (default) or \code{"classic"}.
#' @param alternative Character string specifying the alternative hypothesis
#'   for the Wald tests. One of \code{"two.sided"} (default),
#'   \code{"less"}, or \code{"greater"}.
#' @param ... Further arguments (currently ignored).
#'
#' @details
#' The reported inference is based on asymptotic normality of the coefficient
#' estimates. Regression coefficients are defined as
#' \deqn{\beta = \lambda / (1 - \nu),}
#' where \eqn{\lambda} denotes the dual parameters associated with the
#' non-baseline class and \eqn{\nu} is treated as fixed.
#'
#' Standard errors are obtained from the diagonal of the estimated covariance
#' matrix, and Wald z-statistics are computed as
#' \deqn{z_j = \hat\beta_j / \mathrm{se}(\hat\beta_j).}
#'
#' @return
#' An object of class \code{"summary.GCElogit"}, which is a list with components:
#' \describe{
#'   \item{coefficients}{A matrix with columns \code{Estimate},
#'     \code{Std. Error}, \code{z value}, and \code{Pr(>|z|)}.}
#'   \item{vcov_method}{The covariance estimator used.}
#'   \item{alternative}{The alternative hypothesis used for the Wald tests.}
#' }
#'
#' @seealso
#' \code{\link{confint.GCElogit}},
#'
#' @examples
#' \dontrun{
#' set.seed(7654)
#' betas <- c(-0.5, 1, -0.8, 0.6, 0, 0)
#' sim1 <- DGPlogit(N=1000, beta=betas,
#'                  corr="indep", margins="gaussian")
#' fit <- GCElogit.fit(y=sim1$y, X=sim1$X)
#' summary(fit)
#' }
#' 
#' @importFrom stats printCoefmat
#' @method summary GCElogit
#' @export
summary.GCElogit <- function(object,
                             method = c("sandwich", "classic"),
                             alternative = c("two.sided", "less", "greater"),
                             ...) {
  
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  
  # Covariance estimator
  vc <- switch(
    method,
    sandwich = vcov_GCElogit_sandwich(object),
    classic  = vcov_GCElogit_classic(object)
  )
  
  # Coefficients: beta = lambda / (1 - nu)
  nu   <- as.numeric(object$nu)
  lam  <- as.vector(object$lambda[, 2, drop = FALSE])  # J = 2
  coef <- lam / (1 - nu)
  
  # Standard errors (guard against negative/NA variances)
  vdiag <- diag(vc)
  se <- sqrt(pmax(vdiag, 0))
  
  # Wald statistics
  z <- coef / se
  
  p <- switch(
    alternative,
    two.sided = 2 * pnorm(-abs(z)),
    less      = pnorm(z),
    greater   = pnorm(-z)
  )
  
  # Parameter names
  rn <- rownames(object$lambda)
  if (is.null(rn)) rn <- paste0("x", seq_along(coef))
  nm <- paste0(rn)
  
  # Aliasing (mimic summary.glm: NA for non-estimable terms)
  aliased <- is.na(coef) | is.na(se)
  
  coef_table <- cbind(
    Estimate    = coef,
    `Std. Error` = se,
    `z value`   = z,
    `Pr(>|z|)`  = p
  )
  rownames(coef_table) <- nm
  
  # Assemble output in a glm-like structure
  ans <- list(
    call        = if (!is.null(object$call)) object$call else NULL,
    coefficients = coef_table,
    aliased     = aliased,
    vcov_method = method,
    alternative = alternative
  )
  
  # Keep commonly expected components if available
  if (!is.null(object$terms)) ans$terms <- object$terms
  if (!is.null(object$df.residual)) ans$df.residual <- object$df.residual
  if (!is.null(object$df.null)) ans$df.null <- object$df.null
  if (!is.null(object$deviance)) ans$deviance <- object$deviance
  if (!is.null(object$null.deviance)) ans$null.deviance <- object$null.deviance
  if (!is.null(object$aic)) ans$aic <- object$aic
  
  class(ans) <- "summary.GCElogit"
  ans
}


#' Print Method for \code{summary.GCElogit}
#'
#' Prints a \code{summary.GCElogit} object in a \code{summary.glm}-like format.
#'
#' @param x An object of class \code{"summary.GCElogit"}.
#' @param digits Number of significant digits to use when printing.
#' @param ... Further arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @return The input object \code{x}, invisibly.
#'
#' @method print summary.GCElogit
#' @export
print.summary.GCElogit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  
  cat("\nCall:\n")
  if (!is.null(x$call)) {
    print(x$call)
  } else {
    cat("  <call not available>\n")
  }
  
  cat("\nCoefficients:\n")
  
  ct <- x$coefficients
  
  # Mimic summary.glm formatting
  printCoefmat(ct, digits = digits, signif.stars = getOption("show.signif.stars"), ...)
  
  cat("\nCovariance estimator: ", x$vcov_method, "\n", sep = "")
  
  invisible(x)
}

