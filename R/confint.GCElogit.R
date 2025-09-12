#' Wald confidence intervals for beta = lambda/(1 - nu) from a GCElogit fit
#' S3 method for class "GCElogit"
#'
#' @param object an object of class "gce_mult" returned by gce_mult(), expected to contain:
#'        $lambda (K x J), $hess (Hessian for free params vec(lambda[,2:J])),
#'        and $nu (scalar in (0,1))
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered
#' @param level confidence level (default 0.95)
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
confint.GCElogit <- function(object, parm = NULL, level = 0.95, ...) {
  Vbeta_free <- vcov_GCElogit_sandwich(object)
  nu   <- object$nu
  lamf <- as.vector(object$lambda[, 2, drop=FALSE])  # J=2 -> una sola colonna libera
  betf <- lamf / (1 - nu)
  se   <- sqrt(pmax(diag(Vbeta_free), 0))
  
  # nomi parametri
  rn <- rownames(object$lambda)
  if (is.null(rn)) rn <- paste0("x", seq_along(betf))
  param_names <- paste0("beta[", rn, "]")
  
  # selezione parametri come in stats::confint
  if (is.null(parm)) {
    sel <- seq_along(betf)
  } else if (is.numeric(parm)) {
    sel <- parm
    if (any(sel < 1 | sel > length(betf)))
      stop("parm indices out of range.")
  } else if (is.character(parm)) {
    sel <- match(parm, param_names)
    if (anyNA(sel)) stop("Some 'parm' names not found in parameters.")
  } else {
    stop("parm must be NULL, numeric indices, or character names.")
  }
  
  # calcolo IC
  z <- qnorm(0.5 + level/2)
  lower <- betf[sel] - z * se[sel]
  upper <- betf[sel] + z * se[sel]
  
  out <- cbind(lower = lower, upper = upper)
  rownames(out) <- param_names[sel]
  out
}



# Calcola vcov sandwich per un oggetto "gce_mult"
# Var-cov sandwich per GME logit binario (J=2)
#' @noRd
vcov_GCElogit_sandwich <- function(object, use_sum = TRUE) {
  # Richiesti nell'oggetto:
  #  lambda (K x 2), hess (K x K) sui liberi, nu (scalare),
  #  X (N x K), Y (N) o (N x 2), p (N x 2), e (N x 2)
  stopifnot(!is.null(object$lambda), !is.null(object$hess),
            !is.null(object$nu), !is.null(object$X),
            !is.null(object$Y), !is.null(object$p), !is.null(object$e))
  
  lambda <- object$lambda
  H  <- 0.5 * (object$hess + t(object$hess))   # simmetrizza
  nu <- object$nu
  X  <- object$X                                # N x K
  N  <- nrow(X); K <- ncol(X)
  
  # y_i (successo), pi_i, e_i per la classe non-base (colonna 2)
  if (is.matrix(object$Y)) {
    y <- as.numeric(object$Y[, 2])
  } else {
    y <- as.numeric(object$Y)
  }
  pi <- as.numeric(object$p[, 2])
  e  <- as.numeric(object$e[, 2])
  
  # score per osservazione: S (N x K), riga i = x_i * resid_i
  resid <- -y + pi + e                         # N-vector
  S <- X * resid                               # moltiplicazione per riga (riciclo vettoriale)
  
  # Bread = H^{-1}
  Hinv <- tryCatch(chol2inv(chol(H)), error = function(e) {
    eg <- eigen(H, symmetric = TRUE); d <- pmax(eg$values, 1e-10)
    eg$vectors %*% diag(1/d, length(d)) %*% t(eg$vectors)
  })
  
  # Meat coerente con H: usa SOMMA (default) o MEDIA per entrambi
  if (use_sum) {
    Meat <- crossprod(S)            # = t(S) %*% S = sum_i s_i s_i'
    Vlam <- Hinv %*% Meat %*% Hinv
  } else {
    Meat <- crossprod(S) / N
    Hbar <- H / N
    Hbar_inv <- tryCatch(chol2inv(chol(Hbar)), error = function(e) {
      eg <- eigen(Hbar, symmetric = TRUE); d <- pmax(eg$values, 1e-10)
      eg$vectors %*% diag(1/d, length(d)) %*% t(eg$vectors)
    })
    Vlam <- Hbar_inv %*% Meat %*% Hbar_inv
  }
  
  # Var(beta) con beta = lambda/(1-nu)
  Vbeta <- Vlam / (1 - nu)^2
  Vbeta
}








