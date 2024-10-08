#' Generalized Cross Entropy Linear Regression Models
#'
#' @description Fitting generalized cross entropy (GCE) linear models
#' @author Marco Sandri, Enrico Ciavolino, Maurizio Carpita (\email{gcemodels@gmail.com})
#' @param y numeric, n (Nx1) vector representing the dependent variable where N is the number of observations.
#' @param X numeric, n (NxK) matrix representing a set of independent variables where K is the number of regressors.
#' @param Z numeric, an (KxM) matrix representing support spaces for the regression coefficients (including intercept) where M is the dimension of the support spaces.
#' @param v numeric, an optional argument representing a support space for error terms: 
#' \describe{
#'  \item{(a)}{if missing then \code{v} is a (5x1) vector of equally spaced points in [a,b] interval;}
#'  \item{(b)}{if a scalar (e.g. H) then \code{v} is a (Hx1) vector of equally spaced points in [a,b] interval;}
#'  \item{(c)}{can be a user-supplied vector;}
#'  \item{(d)}{can be a user-supplied matrix.}
#' }
#' Please note that in case (a) and (b) the [a,b] interval is centered around zero, and a and b are calculated using the empirical three-sigma rule Pukelsheim (1994).
#' @param nu numeric, an optional weight parameter representing the trade-off between prediction and precision.
#' @param p0 numeric, optional prior probabilities associated with the regression coefficients.
#' @param w0 numeric, optional prior probabilities associated with the error terms.
#' @param k.sigma numeric, coefficient \code{k} in the k-sigma rule (default \code{k=3}).
#' @param control list, a list of parameters for controlling the fitting process; for \code{GCElm.fit} this is passed to \code{\link{GCElm.control}}.
#' @details Mettere qui eventuali details.
#' @return A \code{list} with the following elements:
#' \itemize{
#'  \item{\code{lambda}, estimated lagrange multipliers;}
#'  \item{\code{coefficients}, regression coefficients;}
#'  \item{\code{var_beta}, variance-covariance matrix of the regression coefficients;}
#'  \item{\code{p}, estimated probabilities associated with the regressions coefficients;}
#'  \item{\code{w}, estimated probabilities associated with the error terms;}
#'  \item{\code{e}, estimated residuals;}
#'  \item{\code{Sp}, the (signal) information of the whole system;}
#'  \item{\code{Sp_k}, the (signal) information associated with the k-th regression coefficient;}
#'  \item{\code{H_p_w}, value of the joint entropies of p and w at the final iteration;}
#'  \item{\code{dH}, delta-H from the Entropy Concentration Theorem;}
#'  \item{\code{ER}, entropy-ratio statistic;}
#'  \item{\code{Pseudo-R2}, pseudo R-squared;}
#'  \item{\code{converged}, convergence (same as in the \code{\link[lbfgs]{lbfgs}} function).}
#' }
#' @references Golan (1996)
#' @examples
#' set.seed(1234)
#' N <- 1000
#' K <- 10
#' betas <- c(-1,1,-1,1,0,0,0,0,0,0)
#' X <- matrix(runif(N*K), nrow = N, ncol = K)
#' y <- X %*% betas + rnorm(N)
#' X <- cbind(rep(1, N), X)
#' Z <- matrix(rep(seq(-10,10,2.5), K+1), nrow = K+1, byrow = TRUE)
#' GCEfit <- GCElm.fit(y, X, Z)
#' coef(GCEfit)      
#' @export
#' @importFrom lbfgs lbfgs
#' @importFrom stats sd
#' @importFrom stats qchisq
#' @importFrom Rcpp evalCpp
#' @useDynLib GCEmodels, .registration=TRUE

GCElm.fit <- function(y, X, Z, v, nu, p0, w0, k.sigma=3, control=GCElm.control()) {
  control <- do.call("GCElm.control", control)
  
  # Check matrix dimensions
  if (ncol(X)!=nrow(Z)) {
    stop(gettextf("matrix Z has %d rows which is not equal to the number of regressors (%d, including intercept)", 
                  NROW(Z), NCOL(X)), domain = NA)
  }
  
  # Check if a column of all ones is present
  if (!any(apply(X,2,function(x) all(x==1)))) {
    X <- cbind(rep(1, nrow(X)), X)
    warning("a column of all ones was added to matrix X (for intercept estimation)")  
  }
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
                     m=control$m, gtol=control$gtol, 
                     max_linesearch=control$max_linesearch,
                     invisible=control$invisible,
                     linesearch_algorithm = control$linesearch_algorithm)
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
  
  coef <- as.vector(beta_hat)
  names(coef) <- dimnames(X)[[2L]]
  info_estim_all <- list(lambda = lambda_hat, coefficients = coef, 
                         var_beta = var_beta, p = p, w = w, e = e, Sp = Sp, S_pk = S_pk, 
                         H_p_w = gce_optim$H, dH = dH, ER = ER, Pseudo_R2 = R2, 
                         converged = gce_optim$convergence, class="GCElm", optim=gce_optim)
  return(info_estim_all)
}

