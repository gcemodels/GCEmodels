#' Linear Regression Model
#'
#' @author Marco Sandri, Enrico Ciavolino, Maurizio Carpita (\email{gcemodels@gmail.com})
#' @param y numeric, n (Nx1) vector representing the dependent variable where N is the number of observations.
#' @param X numeric, n (NxK) matrix representing a set of independent variables where K is number of regressors.
#' @param Z numeric, An (KxM) matrix representing support spaces the for regression coefficients where M is the dimension of the support spaces.
#' @param v An optional argument representing a support space for error terms: (a) if missing then v is a (5x1) vector of equally spaced points in [a,b] interval; (b) if a scalar (e.g. H) then v is a (Hx1) vector of equally spaced points in [a,b] interval; (c) can be a user-supplied vector; (d) can be a user-supplied matrix. Please note that in case (a) and (b) the [a,b] interval is centered around zero, and a and b are calculated using the empirical three-sigma rule Pukelsheim (1994).
#' @param k.sigma Implement the k-sigma rule
#' @param nu numeric, optional: A weight parameter representing the trade-off between prediction and precision.
#' @param p0 numeric, optional: Prior probabilities associated with the regression coefficients
#' @param w0 numeric, optional: Prior probabilities associated with the error terms
#' @param m m
#' @param gtol gtol
#' @param max_linesearch The maximum number of trials for the line search.This parameter controls the number of function and gradients evaluations per iteration for the line search routine. The default value is 20.
#' @param invisible invisible
#' @param linesearch linesearch
#' @details Mettere qui eventuali details.
#' @return A \code{list} with the following elements:
#' @return * \code{lambda}, estimated lagrange multipliers
#' @return * \code{beta}, regression coefficients
#' @return * \code{var_beta}, variance-covariance matrix of the regression coefficients
#' @return * \code{p}, estimated probabilities associated with the regressions coefficients
#' @return * \code{w}, estimated probabilities associated with the error terms
#' @return * \code{e}, estimated residuals
#' @return * \code{Sp}, the (signal) information of the whole system
#' @return * \code{Sp_k}, the (signal) information associated with the k-th regression coefficient
#' @return * \code{H_p_w}, value of the joint entropies of p and w at the final iteration
#' @return * \code{dH}, delta-H from the Entropy Concentration Theorem
#' @return * \code{ER}, entropy-ratio statistic
#' @return * \code{Pseudo-R2}, pseudo R-squared
#' @return * \code{conv}, convergence (same as in the \code{lbfgs} function) 
#' @references Golan (1996)
#' @examples
#' set.seed(1234)
#' N <- 25000
#' K <- 10
#' y <- runif(N)
#' X <- matrix(runif(N*K), nrow = N, ncol = K)
#' X <- cbind(rep(1, N), X)
#' Z <- matrix(rep(c(-1, -0.5, 0, 0.5, 1), K+1), nrow = K+1, byrow = TRUE)
#' GCEfit <- GCElm.fit(y, X, Z, linesearch="LBFGS_LINESEARCH_BACKTRACKING")
#' data.frame(beta = GCEfit$beta,
#'            beta_lb = GCEfit$beta-1.96*sqrt(diag(GCEfit$var_beta)),
#'            beta_ub = GCEfit$beta+1.96*sqrt(diag(GCEfit$var_beta))
#'            )
#' @export
#' @importFrom lbfgs lbfgs
#' @importFrom stats sd
#' @importFrom stats qchisq
#' @importFrom Rcpp evalCpp
#' @useDynLib GCEmodels, .registration=TRUE

GCElm.fit <- function (y, X, Z, v, nu, p0, w0, k.sigma=3, m=6, gtol=0.9, 
                            max_linesearch = 20,
                            invisible=1,
                            linesearch = "LBFGS_LINESEARCH_DEFAULT") {
  dimX <- dim(X)
  N <- dimX[1]
  K <- dimX[2]
  M <- dim(Z)[2]
  if (missing(v)) {
    dimV <- 5
    sd_y <- stats::sd(y)
    v <- seq(from = -k.sigma * sd_y, to = k.sigma * sd_y, length.out = dimV)
  }
  else {
    if (is.vector(v)) {
      len_v <- length(v)
      if (len_v == 1) {
        dimV <- v
        sd_y <- sd(y)
        v <- seq(from = -k.sigma * sd_y, to = k.sigma * sd_y, length.out = dimV)
      }
      else if (len_v > 1) {
        dimV <- len_v
      }
    }
    else if (is.matrix(v) && dim(v)[1] == N) {
      dimV <- dim(v)[2]
    }
  }
  if (missing(nu)) 
    nu <- 0.5
  if (missing(p0)) 
    p0 <- matrix(1/M, nrow = K, ncol = M)
  if (missing(w0)) 
    w0 <- matrix(1/dimV, nrow = N, ncol = dimV)
  env <- new.env()
  env[["X"]] <- X
  env[["y"]] <- y
  env[["Z"]] <- Z
  env[["v"]] <- as.matrix(v)
  env[["nu"]] <- nu
  env[["p0"]] <- p0
  env[["w0"]] <- w0
  env[["N"]] <- N
  env[["K"]] <- K
  env[["M"]] <- M
  env[["J"]] <- dimV		
  lambda0 <- rep(0, N)
  gce_optim <- lbfgs::lbfgs(call_eval=GCElin_objFunct(), 
                     call_grad=GCElin_gradFunct(), 
                     vars = lambda0, environment=env,
                     m=m, gtol=gtol, max_linesearch=max_linesearch,
                     invisible=invisible,
                     linesearch_algorithm = linesearch)
  lambda_hat <- gce_optim$par
  p <- matrix(0, K, M)
  Omega <- rep(0, K)
  for (k in 1:K) {
    temp <- sum(lambda_hat * X[, k])
    for (m in 1:M) {
      p[k, m] <- p0[k, m] * exp(Z[k, m] * temp/(1 - nu))
    }
    Omega[k] <- sum(p[k, ])
    p[k, ] <- p[k, ]/Omega[k]
  }
  beta_hat <- matrix(apply(Z * p, 1, sum), ncol = 1)
  Psi <- rep(0, N)
  if (is.vector(v)) {
    J <- length(v)
    w <- matrix(0, N, J)
    for (n in 1:N) {
      for (j in 1:J) {
        w[n, j] <- w0[n, j] * exp(v[j] * lambda_hat[n]/nu)
      }
      Psi[n] <- sum(w[n, ])
      w[n, ] <- w[n, ]/Psi[n]
    }
    e <- w %*% matrix(v, ncol = 1)
  }
  else if (is.matrix(v) && dim(v)[1] == N) {
    J <- dim(v)[2]
    w <- matrix(0, N, J)
    for (n in 1:N) {
      for (j in 1:J) {
        w[n, j] <- w0[n, j] * exp(v[n, j] * lambda_hat[n]/nu)
      }
      Psi[n] <- sum(w[n, ])
      w[n, ] <- w[n, ]/Psi[n]
    }
    e <- apply(w * v, 1, sum)
  }
  Sp <- sum(p * log(p))/sum(p0 * log(p0))
  S_pk <- rep(0, K)
  for (k in 1:K) {
    S_pk[k] <- sum(p[k, ] * log(p[k, ]))/sum(p0[k, ] * log(p0[k, 
    ]))
  }
  sigma2_beta <- sum(lambda_hat * lambda_hat)/N
  sigma2_e <- rep(0, N)
  if (is.vector(v)) {
    for (n in 1:N) {
      sigma2_e[n] <- sum((v * v) * w[n, ]) - (sum(v * w[n, 
      ]))^2
    }
  }
  else if (is.matrix(v) && dim(v)[1] == N) {
    for (n in 1:N) {
      sigma2_e[n] <- sum((v[n, ] * v[n, ]) * w[n, ]) - 
        (sum(v[n, ] * w[n, ]))^2
    }
  }
  w2_beta <- (sum(1/sigma2_e)/N)^2
  var_beta <- (sigma2_beta/w2_beta) * solve(t(X) %*% X)
  ER <- rep(0, K)
  p_temp <- rep(1/M, M)
  H_R <- -sum(p_temp * log(p_temp))
  for (k in 1:K) {
    H_U <- -sum(p[k, ] * log(p[k, ]))
    ER[k] <- 2 * H_R - 2 * H_U
  }
  R2 <- 1 - Sp
  CC <- K * (M - 1) + N * (J - 2)
  dH <- stats::qchisq(c(0.9, 0.95, 0.99), df = CC)/(2 * N)
  info_estim_all <- list(lambda = lambda_hat, beta = beta_hat, 
                         var_beta = var_beta, p = p, w = w, e = e, Sp = Sp, S_pk = S_pk, 
                         H_p_w = gce_optim$H, dH = dH, ER = ER, Pseudo_R2 = R2, 
                         conv = gce_optim$convergence, optim=gce_optim)
  return(info_estim_all)
}

