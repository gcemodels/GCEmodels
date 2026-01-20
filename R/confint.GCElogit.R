#' Wald confidence intervals for beta = lambda/(1 - nu) from a GCElogit fit
#' S3 method for class "GCElogit"
#'
#' @param object an object of class "gce_mult" returned by gce_mult(), expected to contain:
#'        $lambda (K x J), $hess (Hessian for free params vec(lambda[,2:J])),
#'        and $nu (scalar in (0,1))
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered
#' @param level confidence level (default 0.95)
#' @param method character string specifying the type of variance–covariance estimator to be used. Possible values are \code{"sandwich"} for the robust sandwich estimator based on the score contributions (default), and \code{"classic"} for the classical estimator based on the inverse of the observed Hessian matrix. Any other value will result in an error.
#' @param ... unused, for S3 compatibility
#' @return a matrix with two columns ("lower","upper") and one row per requested parameter
#' @examples
#' set.seed(7654)
#' betas <- c(-0.5, 1, -0.8, 0.6, 0, 0)
#' sim1 <- DGPlogit(N=1000, beta=betas,
#'                  corr="indep", margins="gaussian")
#' fit <- GCElogit.fit(y=sim1$y, X=sim1$X)
#' cbind(betas, fit$coefficients, confint(fit))
#' @importFrom stats qnorm confint
#' @method confint GCElogit
#' @export
#'
#'
confint.GCElogit <- function(object, parm = NULL, level = 0.95,
                             method = "sandwich", ...) {
  # Select covariance estimator
  Vbeta <- switch(
    method,
    sandwich = vcov_GCElogit_sandwich(object),
    classic  = vcov_GCElogit_classic(object),
    stop("Argument 'method' must be either 'sandwich' or 'classic'.")
  )
  
  # Point estimates: beta = lambda / (1 - nu)
  nu    <- as.numeric(object$nu)
  lam   <- as.vector(object$lambda[, 2, drop = FALSE])  # J = 2: single free column
  beta  <- lam / (1 - nu)
  
  # Standard errors
  se <- sqrt(pmax(diag(Vbeta), 0))
  
  # Parameter names
  rn <- rownames(object$lambda)
  if (is.null(rn)) rn <- paste0("x", seq_along(beta))
  param_names <- paste0("beta[", rn, "]")
  
  # Parameter selection (stats::confint semantics)
  if (is.null(parm)) {
    sel <- seq_along(beta)
  } else if (is.numeric(parm)) {
    if (any(parm < 1 | parm > length(beta)))
      stop("Parameter indices out of range.")
    sel <- parm
  } else if (is.character(parm)) {
    sel <- match(parm, param_names)
    if (anyNA(sel))
      stop("Some parameter names were not found.")
  } else {
    stop("Argument 'parm' must be NULL, numeric, or character.")
  }
  # Wald confidence intervals
  z <- qnorm(0.5 + level / 2)
  ci <- cbind(
    lower = beta[sel] - z * se[sel],
    upper = beta[sel] + z * se[sel]
  )
  rownames(ci) <- param_names[sel]

  return(ci)
}


# Classical (model-based) vcov for GCE/GME binary logit (J = 2)
# Returns Var(beta) where beta = lambda / (1 - nu), treating nu as fixed.
#' @noRd
vcov_GCElogit_classic <- function(object, use_sum = TRUE) {
  # Expected fields:
  #   hess : (K x K) Hessian w.r.t. lambda
  #   nu   : scalar
  #   X    : (N x K), used for dimension names
  stopifnot(!is.null(object$hess), !is.null(object$nu), !is.null(object$X))
  
  X  <- object$X
  N  <- nrow(X)
  nu <- as.numeric(object$nu)
  
  # Symmetrise for numerical stability
  H <- 0.5 * (object$hess + t(object$hess))
  
  # Match scaling conventions:
  # - use_sum = TRUE  -> H is a sum over observations
  # - use_sum = FALSE -> H is an average (divide by N)
  H_use <- if (use_sum) H else H / N
  
  # Var(lambda) approximation using a robust inverse
  Hinv <- inv_sym_pd(H_use)
  
  # Delta method for beta = lambda / (1 - nu), with nu treated as constant
  Vbeta <- Hinv / (1 - nu)^2
  dimnames(Vbeta) <- list(colnames(X), colnames(X))
  
  # Enforce numerical symmetry of the covariance matrix  
  asym <- max(abs(Vbeta - t(Vbeta)))
  if (asym > 1e-8) warning("Covariance matrix of regression coefficients shows non-negligible asymmetry (max abs diff = ", asym, ").")
  Vbeta <- 0.5 * (Vbeta + t(Vbeta))

  return(Vbeta)
}


# Sandwich vcov for GCE/GME binary logit (J = 2)
# Returns Var(beta) where beta = lambda / (1 - nu), treating nu as fixed.
#' @noRd
vcov_GCElogit_sandwich <- function(object, use_sum = TRUE, asym_tol = 1e-8) {
  # Expected fields:
  #   hess : (K x K) Hessian w.r.t. lambda
  #   nu   : scalar
  #   X    : (N x K) design matrix
  #   Y    : (N) or (N x 2) response (success indicator in column 2 if matrix)
  #   p    : (N x 2) fitted probabilities (class 2 in column 2)
  #   e    : (N x 2) additional term (class 2 in column 2)
  stopifnot(!is.null(object$hess), !is.null(object$nu), !is.null(object$X),
            !is.null(object$Y),    !is.null(object$p),  !is.null(object$e))
  
  X  <- object$X
  N  <- nrow(X)
  nu <- as.numeric(object$nu)
  
  # Symmetrise for numerical stability
  H <- 0.5 * (object$hess + t(object$hess))
  
  # Extract class-2 quantities
  y  <- if (is.matrix(object$Y)) as.numeric(object$Y[, 2]) else as.numeric(object$Y)
  pi <- as.numeric(object$p[, 2])
  e2 <- as.numeric(object$e[, 2])
  
  # Per-observation score contributions: row i is x_i * (pi_i + e_i - y_i)
  resid <- -y + pi + e2
  S     <- X * resid
  
  # Match scaling conventions: use sums (default) or averages for both bread and meat
  H_use <- if (use_sum) H else H / N
  Meat  <- if (use_sum) crossprod(S) else crossprod(S) / N
  
  # Sandwich covariance for lambda
  Hinv <- inv_sym_pd(H_use)
  Vlam <- Hinv %*% Meat %*% Hinv
  
  # Delta method for beta = lambda / (1 - nu), with nu treated as constant
  Vbeta <- Vlam / (1 - nu)^2
  
  # Enforce numerical symmetry of the covariance matrix of regression coefficients
  asym <- max(abs(Vbeta - t(Vbeta)))
  if (is.finite(asym) && asym > asym_tol) {
    warning("Covariance matrix of regression coefficients shows non-negligible asymmetry (max abs diff = ",
            signif(asym, 6), ").")
  }
  Vbeta <- 0.5 * (Vbeta + t(Vbeta))
  dimnames(Vbeta) <- list(colnames(X), colnames(X))
  
  return(Vbeta)
}


# Invert a symmetric positive-(semi)definite matrix.
# Tries a Cholesky-based inverse; falls back to an eigen-based regularised inverse.
#' @noRd
inv_sym_pd <- function(M, ridge = 1e-10) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M), is.numeric(M), 
            length(ridge) == 1L, ridge > 0)
  # Symmetrise to remove numerical asymmetry
  M <- 0.5 * (M + t(M))
  tryCatch(
    chol2inv(chol(M)),
    error = function(e) {
      eg <- eigen(M, symmetric = TRUE)
      # Regularise small/negative eigenvalues to stabilise inversion
      d_inv <- 1 / pmax(eg$values, ridge)
      # V diag(d_inv) V'
      eg$vectors %*% (d_inv * t(eg$vectors))
    }
  )
}

