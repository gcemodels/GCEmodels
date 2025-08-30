#' @title Generalized Cross Entropy Linear Regression Models
#'
#' @description Fitting generalized cross entropy (GCE) linear models
#' @name GCElm 
#' @author Marco Sandri, Enrico Ciavolino, Maurizio Carpita (\email{gcemodels@gmail.com})
#' @param formula an object of class "formula" (or one that can be coerced to that class); a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which GCElm is called.
#' @param Z numeric, an (KxM) matrix representing support spaces for the regression coefficients (including intercept)  where M is the dimension of the support spaces.
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
#' @param weights an optional vector of ‘prior weights’ to be used in the fitting process; should be NULL or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs; the default is set by the na.action setting of options, and is na.fail if that is unset; the ‘factory-fresh’ default is na.omit; another possible value is NULL, no action; value na.exclude can be useful.
#' @param control list, a list of parameters for controlling the fitting process; for \code{GCElm.fit} this is passed to \code{\link{GCElm.control}}.
#' @param model a logical value indicating whether model frame should be included as a component of the returned value. 
#' @param method the method to be used in fitting the model; the default method \code{GCElm.fit} uses Limited-memory BFGS (L-BFGS); the alternative "model.frame" returns the model frame and does no fitting.
#' @param x logical values indicating whether the model matrix used in the fitting process should be returned as components of the returned value.
#' @param y logical values indicating whether the response vector used in the fitting process should be returned as components of the returned value.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting; this should be NULL or a numeric vector of length equal to the number of cases; one or more \code{\link[stats]{offset}} terms can be included in the formula instead or as well, and if more than one is specified their sum is used; see \code{\link[stats]{model.offset}}.
#' @param contrasts an optional list; see the \code{contrasts.arg} of \code{\link[stats]{model.matrix.default}}.
#' @param ... for GCElm: arguments to be used to form the default control argument if it is not supplied directly; for weights: further arguments passed to or from other methods.
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
#' @return * \code{converged}, convergence (same as in the \code{lbfgs} function) 
#' @references Golan (1996)
#' @examples
#' set.seed(1234)
#' N <- 1000
#' K <- 10
#' betas <- c(-1,1,-1,1,0,0,0,0,0,0)
#' X <- matrix(runif(N*K), nrow = N, ncol = K)
#' y <- X %*% betas + rnorm(N)
#' Z <- matrix(rep(seq(-10,10,2.5), K+1), nrow = K+1, byrow = TRUE)
#' df <- data.frame(y, X)
#' GCEfit <- GCElm(y~., data=df, Z=Z)
#' coef(GCEfit)  
#' @export
#' @importFrom stats model.response
#' @importFrom stats is.empty.model
#' @importFrom stats model.matrix
#' @importFrom stats model.weights
#' @importFrom stats model.offset
#' @importFrom stats .getXlevels
#' @useDynLib GCEmodels, .registration=TRUE 
#' @rawNamespace exportPattern("^[[:alpha:]]+")
#' 
GCElm <- function(formula, data, Z, v, nu, p0, w0, k.sigma=3, weights, subset, na.action, 
                  control = list(), model = TRUE, method = "GCElm.fit", x = FALSE, y = TRUE, 
                  offset, contrasts = NULL, ...){
  cal <- match.call()
  if (missing(data)) 
    data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame")) 
    return(mf)
  if (!is.character(method) && !is.function(method)) 
    stop("invalid 'method' argument")
  if (identical(method, "GCElm.fit")) 
    control <- do.call("GCElm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
  }

  # Default values for v
  dimX <- dim(X)
  N <- dimX[1]
  K <- dimX[2]
  M <- dim(Z)[2]  
  if (missing(v)) {
    dimV <- 5
    sd_y <- stats::sd(Y)
    v <- seq(from = -k.sigma * sd_y, to = k.sigma * sd_y, length.out = dimV)
  }
  else {
    if (is.vector(v)) {
      len_v <- length(v)
      if (len_v == 1) {
        dimV <- v
        sd_y <- sd(Y)
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
  # Default values for nu, p0, and w0  
  if (missing(nu)) 
    nu <- 0.5
  if (missing(p0)) 
    p0 <- matrix(1/M, nrow = K, ncol = M)
  if (missing(w0)) 
    w0 <- matrix(1/dimV, nrow = N, ncol = dimV)    
  
  # Fit the GCE model
  fit <- eval(call(if (is.function(method)) "method" else method, 
                   X=X, y=Y, Z=Z, v=v, nu=nu, p0=p0, w0=w0, k.sigma=k.sigma,
                   control=control))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit <- eval(call(if (is.function(method)) "method" else method, 
                      X=X[, "(Intercept)", drop=FALSE], y=Y, Z=Z, v=v, nu=nu, p0=p0, w0=w0, k.sigma=k.sigma,
                      control=control))
  #fit$null.deviance <- fit2$deviance
  }
  if (fit$converged!=0)
    warning(paste0("The optimization algorithm did not converged (code ", fit$converged,")"))  
  if (model) 
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x) 
    fit$x <- X
  if (!y) 
    fit$y <- NULL
  
  qrX <- base::qr(X)
  rankX <- qrX$rank
  df_resid <- NROW(X) - rankX
  
  Yhat <- fit$fitted.values
  
  structure(c(fit, list(call = cal, formula = formula, terms = mt, 
                        data = data, offset = offset, control = control, method = method, 
                        fitted.values = Yhat, residuals = Yhat - Y,
                        qr = qrX, df.residual = df_resid, rank = rankX, nu = nu, 
                        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf))), 
            class = c(fit$class, c("lm")))
}


