#' Auxiliary for Controlling GCE Fitting
#' 
#' @description Auxiliary function for GCE fitting. Typically only used internally by \code{\link{GCElm.fit}}.
#' @author Marco Sandri, Enrico Ciavolino, Maurizio Carpita (\email{gcemodels@gmail.com})
#' @param m numeric, the number of corrections to approximate the inverse Hessian matrix; see \code{\link[lbfgs]{lbfgs}}.
#' @param gtol numeric, a parameter to control the accuracy of the line search routine; see \code{\link[lbfgs]{lbfgs}}.
#' @param max_linesearch numeric, the maximum number of trials for the line search; this parameter controls the number of function and gradients evaluations per iteration for the line search routine; the default value is 20; see \code{\link[lbfgs]{lbfgs}}.
#' @param invisible numeric, defaults to 0; set to 1 to suppress console output.
#' @param linesearch_algorithm string, the line search algorithm; this parameter specifies a line search algorithm to be used by the L-BFGS routine; see \code{\link[lbfgs]{lbfgs}}.
#' @details The control argument of GCElm.fit is by default passed to the control argument of GCElm.fit, which uses its elements as arguments to GCElm.control: the latter provides defaults and sanity checking.
#' @return A \code{list} with components named as the arguments.
#' @references Golan (1996)
#' @examples
#' set.seed(1234)
#' N <- 25000
#' K <- 10
#' y <- runif(N)
#' X <- matrix(runif(N*K), nrow = N, ncol = K)
#' X <- cbind(rep(1, N), X)
#' Z <- matrix(rep(c(-1, -0.5, 0, 0.5, 1), K+1), nrow = K+1, byrow = TRUE)
#' GCElm_control <- GCElm.control(m=8, gtol=0.95)
#' GCEfit <- GCElm.fit(y, X, Z, control=GCElm_control)
#' data.frame(beta = GCEfit$beta,
#'            beta_lb = GCEfit$beta-1.96*sqrt(diag(GCEfit$var_beta)),
#'            beta_ub = GCEfit$beta+1.96*sqrt(diag(GCEfit$var_beta))
#'            )
#' @export

GCElm.control <- function (m = 6, gtol = 0.9, max_linesearch = 20, invisible = 1, linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING") {
  if (!is.numeric(m) || m <= 0) 
    stop("value of 'm' must be > 0")
  if (!is.numeric(gtol) || gtol <= 0) 
    stop("value of 'gtol' must be > 0")
  if (!is.numeric(max_linesearch) || max_linesearch <= 0) 
    stop("value of 'max_linesearch' must be > 0")  
  if (!is.numeric(invisible) || (invisible!=0 && invisible!=1)) 
    stop("invisible must be 0 or 1")
  list(m=m, gtol=gtol, max_linesearch=max_linesearch, invisible=invisible, linesearch_algorithm=linesearch_algorithm)
}
