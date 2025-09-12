#' @title Fit a binary generalized cross-entropy (GCE) logit model
#'
#' @description
#' \code{GCElogit.fit()} estimates a two-class (J = 2) multinomial-logit–type model
#' under the generalized cross-entropy (GCE) framework. The first column of the
#' coefficient matrix is fixed to zero (baseline), so only the second column is
#' optimized. Class posterior probabilities and latent errors are estimated
#' jointly using a finite error support \code{v}.
#'
#' @param y A binary response vector of length \eqn{N} with values in \code{0/1}.
#'   If a binary \code{factor} is supplied, it is converted to \code{0/1} with a
#'   warning: the first level is coded as \code{0}, the second as \code{1}.
#'   Passing a matrix is not allowed.
#' @param X A numeric design matrix of size \eqn{N \times K}. No column may be
#'   constant. The intercept (column of ones) is \emph{not} added automatically.
#' @param v Error support. If missing, the default is
#'   \code{seq(-1, 1, length.out = 5) / sqrt(N)}. If a single numeric is passed,
#'   it is interpreted as \code{length.out} for the default grid. If a numeric
#'   vector is passed, it is used as-is.
#' @param nu Mixing parameter in \eqn{(0,1)} controlling the trade-off between
#'   the two entropy components; default \code{0.5}.
#' @param p0 Optional prior matrix for class probabilities of size \eqn{N \times 2}.
#'   If missing, uniform priors are used.
#' @param w0 Optional prior array for error weights of size \eqn{N \times 2 \times M},
#'   where \eqn{M = length(v)}. If missing, uniform priors are used.
#' @param optim_control A list collecting optimization options. Exactly three
#'   entries are handled by \code{GCElogit.fit()}:
#'   \itemize{
#'     \item \code{method}: optimization method for \code{\link[stats]{optim}}
#'       (e.g., \code{"BFGS"}, \code{"CG"}, \code{"Nelder-Mead"}, \code{"L-BFGS-B"}).
#'       Default \code{"BFGS"}.
#'     \item \code{use_grad}: logical, whether to use an analytic gradient.
#'       Default \code{TRUE}.
#'     \item \code{grad_fun}: an alternative gradient function with the same
#'       signature as \code{GCElogit_gradFunct}, or \code{NULL} to use
#'       \code{GCElogit_gradFunct}. Default \code{NULL}.
#'   }
#'   All other elements supplied in \code{optim_control} are passed \emph{verbatim}
#'   to the \code{control} argument of \code{\link[stats]{optim}} (e.g.,
#'   \code{maxit}, \code{reltol}, \code{trace}, \code{parscale}, \code{fnscale}, ...).
#'   No additional checks are performed; \code{optim()} handles them.
#'
#' @details
#' The model fixes the first column of \eqn{\Lambda} to zero (identification) and
#' optimizes the second column only. With linear index \eqn{X \Lambda}, posterior
#' class probabilities \eqn{p} are proportional to
#' \deqn{p \propto p_0 \odot \exp\{ (X \Lambda) / (1 - \nu) \},}
#' normalized by rows. For each observation and class, latent error weights over
#' the support \code{v} are updated proportionally to
#' \deqn{w_{ijm} \propto w_{0,ijm} \exp\{ (X\Lambda)_{ij} \, v_m / \nu \},}
#' then normalized over \eqn{m}, yielding error means
#' \eqn{e_{ij} = \sum_m v_m w_{ijm}}.
#'
#' The objective minimized by \code{optim()} is:
#' \deqn{
#' \mathcal{L}(\Lambda) =
#'   -\sum_{k,j}\lambda_{kj} (X^\top Y)_{kj}
#'   + (1-\nu)\sum_i \log \Omega_i
#'   + \nu \sum_{i,j} \log \Psi_{ij},
#' }
#' where \eqn{\Omega_i = \sum_j p_{ij}} and \eqn{\Psi_{ij} = \sum_m w_{ijm}}.
#' The analytic gradient \code{GCElogit_gradFunct()} is used when \code{use_grad = TRUE}.
#'
#' Marginal effects are computed using the general formula
#' \eqn{\partial p_{ij} / \partial x_k = (p_{ij} / (1-\nu))(\lambda_{kj} - \sum_r \lambda_{kr} p_{ir})}
#' and averaged over observations. In the binary case with \eqn{\lambda_{\cdot,1} = 0}
#' and \eqn{\lambda_{\cdot,2} = \beta}, this simplifies to
#' \eqn{\partial p_{i2}/\partial x_k = p_{i1} p_{i2} \beta_k / (1-\nu)} and
#' \eqn{\partial p_{i1}/\partial x_k = -\partial p_{i2}/\partial x_k}.
#'
#' @return A list with class \code{"GCElogit"} containing:
#' \itemize{
#'   \item \code{coefficients}: \eqn{K}-vector of slope coefficients \eqn{\beta = \Lambda_{\cdot,2}/(1-\nu)}.
#'   \item \code{hess}: Hessian matrix returned by \code{optim()}.
#'   \item \code{p}: \eqn{N \times 2} matrix of posterior class probabilities.
#'   \item \code{w}: \eqn{N \times 2 \times M} array of error weights.
#'   \item \code{e}: \eqn{N \times 2} matrix of error means.
#'   \item \code{v}: numeric vector with the error support actually used.
#'   \item \code{nu}: mixing parameter.
#'   \item \code{marg_eff}: \eqn{2 \times K} matrix of average marginal effects (rows = classes).
#'   \item \code{Sp}, \code{S_p_i}, \code{p_e_i}: information measures (entropy-based).
#'   \item \code{H_p_w}: negative optimized objective value (entropy of \eqn{p} and \eqn{w} parts).
#'   \item \code{ER}: entropy ratio statistic.
#'   \item \code{Pseudo_R2}: pseudo-\eqn{R^2} defined as \eqn{1 - Sp}.
#'   \item \code{CM}: \eqn{2 \times 2} confusion matrix based on \code{argmax(p + e)}.
#'   \item \code{optim_convergence}: convergence code from \code{optim()}.
#'   \item \code{X}, \code{y}: data used for fitting.
#' }
#'
#' @section Warnings:
#' If \code{y} is a factor with two levels, it is converted to \code{0/1} with a warning,
#' where \code{levels(y)[1]} is coded as \code{0} and \code{levels(y)[2]} as \code{1}.
#'
#' @references
#' Golan, A., Judge, G., & Miller, D. (1996). \emph{Maximum Entropy Econometrics: Robust Estimation with Limited Data}. Wiley.  
#' Golan, A. (1988). Information and entropy econometrics — A review and synthesis.
#'
#' @seealso \code{\link[stats]{optim}}
#'
#' @examples
#' set.seed(7654)
#' sim1 <- DGPlogit(N=1000, beta=c(-0.5, 1, -0.8, 0.6, 0, 0),
#'                  corr="indep", margins="gaussian")
#' fit <- GCElogit.fit(y=sim1$y, X=sim1$X)
#' fit$coefficients
#' @export
GCElogit.fit <- function(y, X, v,
                         nu = 0.5, p0, w0,
                         optim_control = list(method = "BFGS",
                                              use_grad = TRUE,
                                              grad_fun = NULL)) {
  
  # --- Check X
  if (!is.matrix(X)) stop("X must be a numeric matrix.")
  if (!is.numeric(X)) stop("X must be numeric.")
  if (any(apply(X, 2, var) == 0)) stop("No column of X can be constant.")
  X <- cbind(1, X)
  N <- nrow(X); K <- ncol(X); J <- 2
  
  # --- Check y
  if (is.matrix(y)) stop("y must be a 0/1 vector, not a matrix.")
  if (is.factor(y)) {
    if (nlevels(y) != 2) stop("y must have exactly two levels if it is a factor.")
    lvl <- levels(y)
    warning(sprintf(
      "Factor y converted to numeric 0/1. Category '%s' coded as 0, category '%s' coded as 1.",
      lvl[1], lvl[2]
    ))
    y <- as.numeric(y) - 1
  }
  if (!is.numeric(y)) stop("y must be numeric (0/1) or a binary factor.")
  if (length(y) != N) stop("Length of y does not match the number of rows in X.")
  if (!all(y %in% c(0, 1))) stop("y must contain only 0 and 1 values.")
  Y <- cbind(1 - y, y)
  
  # --- Handle v
  if (missing(v)) {
    v <- seq(from = -1, to = 1, length.out = 5) / sqrt(N)
  } else if (length(v) == 1 && is.numeric(v)) {
    v <- seq(from = -1, to = 1, length.out = as.integer(v)) / sqrt(N)
  } else if (is.numeric(v) && length(v) > 1) {
    v <- as.numeric(v)
  } else {
    stop("v must be either a single numeric value (length of support) or a numeric vector.")
  }
  M <- length(v)
  
  # --- Priors
  p0_is_uniform <- FALSE
  if (missing(p0)) { p0 <- matrix(1/J, N, J); p0_is_uniform <- TRUE }
  if (missing(w0)) { w0 <- array(1/M, dim = c(N, J, M)) }
  
  # --- Gradient
  if (optim_control$use_grad) {
    grad_to_use <- if (is.null(optim_control$grad_fun)) GCElogit_gradFunct else grad_fun
  } else {
    grad_to_use <- NULL
  }
  
  # --- Extract optimization options
  if (is.null(optim_control) || !is.list(optim_control)) {
    stop("optim_control must be a list.")
  }
  # --- Extract special parameters from optim_control
  method   <- if (is.null(optim_control$method))  "BFGS" else optim_control$method
  use_grad <- if (is.null(optim_control$use_grad)) TRUE   else optim_control$use_grad
  grad_fun <- if (is.null(optim_control$grad_fun)) NULL   else optim_control$grad_fun
  
  # --- Everything else goes to control
  special  <- c("method", "use_grad", "grad_fun")
  control  <- optim_control[setdiff(names(optim_control), special)]
  
  # --- Optimization (classic call to optim)
  lambda0 <- rep(0, K)
  gce_optim <- optim(
    par     = lambda0,
    fn      = GCElogit_objFunct,
    gr      = grad_to_use,
    method  = method,
    control = control,
    hessian = TRUE,
    Y = Y, X = X, v = v, nu = nu, p0 = p0, w0 = w0,
    N = N, K = K, J = J, M = M
  )
  
  # --- Results
  lambda <- matrix(gce_optim$par, K, J - 1)
  lambda <- cbind(rep(0, K), lambda)
  
  temp <- X %*% lambda
  p <- p0 * exp(temp / (1 - nu))
  p <- p / rowSums(p)
  
  w <- array(0, dim = c(N, J, M))
  e <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      w[i, j, ] <- w0[i, j, ] * exp(temp[i, j] * v / nu)
      w[i, j, ] <- w[i, j, ] / sum(w[i, j, ])
      e[i, j]   <- sum(v * w[i, j, ])
    }
  }
  
  # Marginal effects
  dp_dx <- array(0, c(N, J, K))
  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:K) {
        dp_dx[i, j, k] <- (p[i, j] / (1 - nu)) * (lambda[k, j] - sum(lambda[k, ] * p[i, ]))
      }
    }
  }
  ave_dp_dx <- apply(dp_dx, c(2, 3), mean)
  
  # Confusion matrix
  y_hat <- p + e
  CM    <- matrix(0, J, J)
  for (n in 1:N) {
    t1 <- match(1, Y[n, ])
    t2 <- which.max(y_hat[n, ])
    CM[t1, t2] <- CM[t1, t2] + 1
  }
  
  # Information measures, Golan (1988)
  if (p0_is_uniform == TRUE) {
    Sp <- -sum(p * log(p)) / (N * log(J))
  } else {
    Sp <-  sum(p * log(p)) / sum(p0 * log(p0))
  }
  
  S_p_i <- rep(0, N)
  for (i in 1:N) {
    if (p0_is_uniform == TRUE) {
      S_p_i[i] <- -sum(p[i, ] * log(p[i, ])) / log(J)
    } else {
      S_p_i[i] <-  sum(p[i, ] * log(p[i, ])) / sum(p0[i, ] * log(p0[i, ]))
    }
  }
  
  # Error bounds based on the Fano's [weaker] inequality
  p_e_i <- rep(0, N)
  for (i in 1:N) {
    p_e_i[i] <- S_p_i[i] - log(J)
  }
  
  # Entropy ratio Statistic
  if (p0_is_uniform == TRUE) {
    ER <- 2 * N * log(J) * (1 - Sp)
  } else {
    ER <- 2 * (-sum(p0[i, ] * log(p0[i, ]))) * (1 - Sp)
  }
  
  # Pseudo R-squared
  R2 <- 1 - Sp
  
  # Coefficients
  betas <- lambda[,2]/(1-nu)
  
  out <- list(
    coefficients        = betas,
    lambda              = lambda,
    hess                = gce_optim$hessian,
    p                   = p,
    w                   = w,
    e                   = e,
    v                   = v,
    nu                  = nu,
    marg_eff            = ave_dp_dx,
    Sp                  = Sp,
    S_p_i               = S_p_i,
    p_e_i               = p_e_i,
    H_p_w               = -gce_optim$value,
    ER                  = ER,
    Pseudo_R2           = R2,
    CM                  = CM,
    optim_convergence   = gce_optim$convergence,
    X                   = X,
    Y                   = Y,
    y                   = y
  )
  
  class(out) <- c("GCElogit", class(out))
  return(out)
}




#' @keywords internal
#' @noRd
GCElogit_objFunct <- function(lambda_vector, Y, X, v, nu, p0, w0, N, K, J, M) {
  
  lambda <- matrix(lambda_vector, K, (J - 1))
  lambda <- cbind(rep(0, K), lambda)
  
  temp  <- X %*% lambda
  p <- p0 * exp(temp / (1 - nu))
  Omega <- apply(p, 1, sum)
  
  w <- array(0, c(N, J, M))
  Psi <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      w[i, j, ] <- w0[i, j, ] * exp(temp[i, j] * v / nu)
      Psi[i, j] <- sum(w[i, j, ])
    }
  }
  
  like  <- - sum(lambda * (t(X) %*% Y)) + (1 - nu) * sum(log(Omega)) + nu * sum(log(Psi))
  return(like)
  
}


#' @keywords internal
#' @noRd
GCElogit_gradFunct <- function(lambda_vector, Y, X, v, nu, p0, w0, N, K, J, M) {
  
  lambda <- matrix(lambda_vector, K, (J - 1))
  lambda <- cbind(rep(0, K), lambda)
  
  temp  <- X %*% lambda
  p <- p0 * exp(temp / (1 - nu))
  Omega <- apply(p, 1, sum)
  for (i in 1:N) {
    p[i, ] <- p[i, ] / Omega[i]
  }
  
  w   <- array(0, c(N, J, M))
  Psi <- matrix(0, N, J)
  e   <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      w[i, j, ] <- w0[i, j, ] * exp(temp[i, j] * v / nu)
      Psi[i, j] <- sum(w[i, j, ])
      w[i, j, ] <- w[i, j, ] / Psi[i, j]
      e[i, j]   <- sum(v * w[i, j, ])
    }
  }
  
  dlambda <- rep(0, K * (J -1) )
  index_lambda <- 1
  for (j in 2:J) {
    for (k in 1:K) {
      dlambda[index_lambda] <- - sum(Y[, j] * X[, k]) +
        sum(p[, j] * X[, k]) +
        sum(e[, j] * X[, k])
      index_lambda <- index_lambda + 1
    }
  }
  return(dlambda)
  
}




