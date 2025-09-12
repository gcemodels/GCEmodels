#' Data Generating Process for Binary Logistic Regression
#'
#' @description
#' \code{DGPlogit()} simulates binary outcomes from a logistic regression model
#' with flexible correlation structures and user-specified marginal distributions
#' for predictors. Predictors are generated via a Gaussian copula to enforce a
#' desired correlation matrix while allowing heterogeneous margins.
#'
#' @param N Integer. Number of observations.
#' @param beta Numeric vector of coefficients including the intercept
#'   (length \eqn{K+1}). The first element is the intercept, followed by
#'   slopes for each of the \eqn{K} predictors.
#' @param corr Correlation structure among predictors. One of:
#'   \code{"indep"} (independent predictors),
#'   \code{"ar1"} (AR(1) with parameter \code{rho}),
#'   \code{"exp"} (exponential with decay \code{phi}),
#'   or \code{"user"} (user-specified correlation matrix \code{R}).
#' @param rho Numeric scalar in (-1, 1). Correlation parameter for AR(1).
#' @param phi Positive scalar. Decay parameter for exponential correlation.
#' @param R User-supplied correlation matrix (\eqn{K \times K}) if
#'   \code{corr = "user"}.
#' @param mu Numeric vector of means for margins (recycled to length K).
#'   For Gaussian/t margins, used as mean; for lognormal, used as \code{meanlog}.
#' @param sd Numeric vector of standard deviations (recycled to length K).
#'   For Gaussian/t margins, used as scale; for lognormal, used as \code{sdlog}.
#' @param margins Character vector of length K specifying marginal distributions
#'   for each predictor. Options: \code{"gaussian"}, \code{"t"},
#'   \code{"lognormal"}, \code{"uniform"}, \code{"bernoulli"}.
#' @param df Degrees of freedom for t margins (scalar or length K).
#' @param unif_min,unif_max Lower and upper bounds for uniform margins
#'   (scalars or length K).
#' @param bernoulli_p Success probability for Bernoulli margins
#'   (scalar or length K).
#' @param standardize Logical. If \code{TRUE}, continuous predictors are centered
#'   and scaled to unit variance.
#' @param include_intercept_col Logical. If \code{TRUE}, an intercept column of 1
#'   is prepended to \code{X}.
#' @param target_prevalence Optional numeric scalar in (0, 1). If provided,
#'   the intercept is shifted so that the simulated prevalence
#'   \eqn{E[y]} matches \code{target_prevalence} (within tolerance).
#' @param intercept_bounds Numeric length-2 vector. Interval in which to search
#'   for the intercept shift when calibrating prevalence.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @details
#' The function generates latent Gaussian vectors \eqn{Z_i \sim N(0, R)} for
#' each observation, where \eqn{R} is the correlation matrix implied by
#' \code{corr}, \code{rho}, \code{phi}, or provided as \code{R}.
#' These are transformed to uniform \eqn{U_i = \Phi(Z_i)} and then mapped into
#' the specified marginal distributions via their inverse CDFs.
#'
#' The linear predictor is
#' \deqn{\eta = \beta_0 + X \beta,}
#' with success probability
#' \deqn{p = \frac{1}{1 + \exp(-\eta)}.}
#' Outcomes \code{y} are drawn as \code{Bernoulli(p)}.
#'
#' If \code{target_prevalence} is provided, the intercept \eqn{\beta_0}
#' is shifted by solving
#' \deqn{\frac{1}{N}\sum_i \text{logit}^{-1}(\eta_i + \delta) = \text{target\_prevalence}}
#' for \eqn{\delta} via root-finding within \code{intercept_bounds}.
#'
#' @return A list of class \code{"DGPlogit"} with components:
#' \itemize{
#'   \item \code{X}: Simulated design matrix (\eqn{N \times K}, or \eqn{N \times (K+1)}
#'     if \code{include_intercept_col=TRUE}).
#'   \item \code{y}: Binary outcome vector of length N.
#'   \item \code{p}: Vector of success probabilities.
#'   \item \code{eta}: Linear predictor.
#'   \item \code{beta}: Coefficient vector (with intercept).
#'   \item \code{slopes}, \code{intercept}, \code{intercept_shift}.
#'   \item \code{R_used}: Correlation matrix actually used.
#'   \item \code{margins}, \code{mu}, \code{sd}, \code{df}, \code{unif_min},
#'     \code{unif_max}, \code{bernoulli_p}: Margin specifications.
#'   \item \code{corr}, \code{rho}, \code{phi}: Correlation info.
#'   \item \code{standardize}, \code{include_intercept_col}, \code{target_prevalence}.
#'   \item \code{seed}: RNG seed if set.
#' }
#'
#' @examples
#' # Example 1: Gaussian independent predictors
#' set.seed(1)
#' sim1 <- DGPlogit(N=1000, beta=c(-0.5, 1, -0.8, 0.6),
#'                  corr="indep", margins="gaussian")
#'
#' # Example 2: AR(1) correlation, mixed margins
#' sim2 <- DGPlogit(N=500, beta=c(0.2, 0.5, -1, 0.8),
#'                  corr="ar1", rho=0.5,
#'                  margins=c("gaussian","t","bernoulli"),
#'                  df=5, bernoulli_p=0.3, standardize=TRUE)
#'
#' # Example 3: Target prevalence calibration
#' sim3 <- DGPlogit(N=1000, beta=c(0, 1, -1),
#'                  corr="exp", phi=0.7,
#'                  target_prevalence=0.25)
#'#' @seealso \code{\link[stats]{optim}}, \code{\link[stats]{rbinom}}, \code{\link[stats]{qnorm}}
#'
#' @export
#' @importFrom stats optim pnorm qlnorm qnorm qunif rbinom rnorm uniroot
#' 
#' 
# DGP for binary logistic regression with flexible covariates and correlation
DGPlogit <- function(
    N,
    beta,                      # coefficients INCLUDING intercept
    corr = c("indep", "ar1", "exp", "user"),
    rho = NULL,                # AR(1) parameter in (-1, 1)
    phi = NULL,                # exponential decay (>0): R_ij = exp(-phi * |i-j|)
    R   = NULL,                # user-provided correlation matrix (K x K)
    mu  = 0,                   # means for margins (recycled to length K)
    sd  = 1,                   # sds for margins (recycled to length K, for continuous)
    margins = "gaussian",      # vector of length K: gaussian/t/lognormal/uniform/bernoulli
    df = 5,                    # df for t margins (scalar or length K)
    unif_min = 0,              # lower bounds for uniform margins (scalar or length K)
    unif_max = 1,              # upper bounds for uniform margins (scalar or length K)
    bernoulli_p = 0.5,         # success prob for Bernoulli margins (scalar or length K)
    standardize = FALSE,       # if TRUE, center and scale continuous X columns
    include_intercept_col = FALSE,  # if TRUE, return X with a leading 1 column
    target_prevalence = NULL,  # if not NULL, adjust intercept to hit E[y]=target
    intercept_bounds = c(-20, 20),  # search bounds for intercept adjustment
    seed = NULL                # RNG seed
) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- Dimensions
  if (length(beta) < 2) stop("beta must include an intercept and at least one slope.")
  K <- length(beta) - 1
  if (K < 1) stop("Need at least one predictor (length(beta) >= 2).")
  
  # --- Arg checks (basic)
  corr <- match.arg(corr)
  # recycle helpers
  recycle_to <- function(x, n) if (length(x) %in% c(1, n)) rep(x, length.out = n) else
    stop(sprintf("Length mismatch: expected length 1 or %d, got %d.", n, length(x)))
  margins   <- recycle_to(tolower(margins), K)
  mu        <- recycle_to(mu, K)
  sd        <- recycle_to(sd, K)
  df        <- recycle_to(df, K)
  unif_min  <- recycle_to(unif_min, K)
  unif_max  <- recycle_to(unif_max, K)
  bernoulli_p <- recycle_to(bernoulli_p, K)
  
  # --- Correlation matrix R_used
  if (corr == "user") {
    if (is.null(R)) stop("corr='user' requires a correlation matrix R.")
    if (!is.matrix(R) || any(dim(R) != c(K, K))) stop("R must be K x K.")
    if (any(abs(R - t(R)) > 1e-8)) stop("R must be symmetric.")
    if (any(diag(R) != 1)) stop("Diagonal of R must be all ones.")
    R_used <- R
  } else if (corr == "ar1") {
    if (is.null(rho)) stop("corr='ar1' requires rho in (-1,1).")
    if (!is.numeric(rho) || length(rho) != 1L || abs(rho) >= 1) stop("rho must be a scalar in (-1,1).")
    idx <- seq_len(K)
    R_used <- outer(idx, idx, function(i, j) rho^abs(i - j))
  } else if (corr == "exp") {
    if (is.null(phi)) stop("corr='exp' requires positive phi.")
    if (!is.numeric(phi) || length(phi) != 1L || phi <= 0) stop("phi must be a positive scalar.")
    idx <- seq_len(K)
    R_used <- outer(idx, idx, function(i, j) exp(-phi * abs(i - j)))
  } else if (corr == "indep") {
    R_used <- diag(K)
  } else {
    stop("Unknown corr type.")
  }
  
  # --- Construct Gaussian copula latent Z ~ N(0, R_used)
  # Cholesky may fail if R_used is near-singular; add tiny jitter if needed
  chol_safe <- function(R) {
    out <- try(chol(R), silent = TRUE)
    if (inherits(out, "try-error")) {
      eps <- 1e-10
      out <- chol(R + eps * diag(nrow(R)))
    }
    out
  }
  L <- chol_safe(R_used)
  Z <- matrix(rnorm(N * K), N, K) %*% L             # N x K latent Gaussian
  U <- pnorm(Z)                                     # Gaussian copula U(0,1)
  
  # --- Map to requested margins (vectorized)
  X <- matrix(NA_real_, N, K)
  for (j in seq_len(K)) {
    mj <- margins[j]
    if (mj == "gaussian") {
      X[, j] <- qnorm(U[, j], mean = mu[j], sd = sd[j])
    } else if (mj == "t") {
      X[, j] <- qt(U[, j], df = df[j]) * sd[j] + mu[j]
    } else if (mj == "lognormal") {
      # interpret mu/sd as meanlog/sdlog for lognormal
      X[, j] <- qlnorm(U[, j], meanlog = mu[j], sdlog = sd[j])
    } else if (mj == "uniform") {
      X[, j] <- qunif(U[, j], min = unif_min[j], max = unif_max[j])
    } else if (mj == "bernoulli") {
      X[, j] <- as.numeric(U[, j] < bernoulli_p[j])
    } else {
      stop(sprintf("Unsupported margin '%s'. Allowed: gaussian, t, lognormal, uniform, bernoulli.", mj))
    }
  }
  
  # --- Optional standardization (continuous columns only)
  if (isTRUE(standardize)) {
    cont <- !(margins %in% "bernoulli")
    if (any(cont)) {
      Xc <- scale(X[, cont, drop = FALSE])
      # keep attributes? remove to return plain numeric matrix
      attr(Xc, "scaled:center") <- NULL
      attr(Xc, "scaled:scale")  <- NULL
      X[, cont] <- Xc
    }
  }
  
  # --- Linear predictor, probability, and Bernoulli outcome
  intercept <- beta[1]
  slopes    <- beta[-1]           # length K
  
  eta <- drop(intercept + X %*% slopes)
  sigmoid <- function(x) 1 / (1 + exp(-x))
  p <- sigmoid(eta)
  
  # --- Optional: calibrate intercept to reach target prevalence
  intercept_shift <- 0
  if (!is.null(target_prevalence)) {
    if (!is.numeric(target_prevalence) || length(target_prevalence) != 1 ||
        target_prevalence <= 0 || target_prevalence >= 1) {
      stop("target_prevalence must be a scalar in (0,1).")
    }
    f <- function(delta) mean(sigmoid(eta + delta)) - target_prevalence
    root <- try(uniroot(f, interval = intercept_bounds), silent = TRUE)
    if (inherits(root, "try-error")) {
      warning("Could not calibrate intercept to target_prevalence within bounds; using original intercept.")
    } else {
      intercept_shift <- root$root
      eta <- eta + intercept_shift
      p   <- sigmoid(eta)
    }
  }
  
  y <- rbinom(N, size = 1, prob = p)
  
  # --- Optionally add an intercept column to X (for convenience)
  X_out <- if (include_intercept_col) cbind(`(Intercept)` = 1, X) else X
  
  # --- Return
  out <- list(
    X                     = X_out,
    y                     = y,
    p                     = p,
    eta                   = eta,
    beta                  = beta,
    slopes                = slopes,
    intercept             = intercept + intercept_shift,
    intercept_shift       = intercept_shift,
    margins               = margins,
    mu                    = mu,
    sd                    = sd,
    df                    = df,
    unif_min              = unif_min,
    unif_max              = unif_max,
    bernoulli_p           = bernoulli_p,
    R_used                = R_used,
    corr                  = corr,
    rho                   = rho,
    phi                   = phi,
    standardize           = standardize,
    include_intercept_col = include_intercept_col,
    target_prevalence     = target_prevalence,
    seed                  = seed
  )
  class(out) <- c("DGPlogit", class(out))
  return(out)
}
