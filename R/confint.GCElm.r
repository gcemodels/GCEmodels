#' Confidence Intervals for Model Parameters
#'
#' Computes confidence intervals for one or more parameters in a fitted
#' model. This is a method for objects of class \code{lm}.
#'
#' @param object an object of class \code{lm}.
#' @param parm a specification of which parameters are to be given
#'        confidence intervals, either a vector of numbers or a vector
#'        of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The intervals are based on the asymptotic normality of the estimates
#' and use the estimated covariance matrix of the coefficients.
#'
#' @return A matrix with columns giving lower and upper confidence limits
#' for each parameter.
#'
#' @seealso
#' \code{\link[stats]{coef}}, \code{\link[GCEmodels]{GCElm}}
#'
#' @examples
#' fit <- lm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris)
#' confint(fit)
#'
#' @importFrom stats qt
#' @importFrom stats setNames
#' @export
confint.GCElm <- function (object, parm, level = 0.95, ...) 
{
  cf <- coef(object)
  # ses <- sqrt(diag(vcov(object)))
  ses <- sqrt(diag(object$var_beta))/object$nu
  pnames <- names(ses)
  if (is.matrix(cf)) 
    cf <- setNames(as.vector(cf), pnames)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qt(a, object$df.residual)
  #pct <- .format_perc(a, 3)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                                   pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}

