#' GCEmodels: Generalized (Cross) Entropy Models for Linear Regression
#'
#' @description
#' Strumenti per la stima e l'analisi di modelli lineari basati su approcci di entropia
#' (GME/GCE), con funzioni di fitting, summary e utilità per diagnostica e confronto.
#'
#' @details
#' Funzioni principali:
#' - [GCElm()] per la stima del modello
#' - [summary.GCElm()] per il riassunto dei risultati
#'
#' @references
#' Judge, G. & Mittelhammer, R. (2012). *An Information Theoretic Approach to Econometrics*. Academic Press.  
#' Golan, A., Judge, G., & Miller, D. (1996). *Maximum Entropy Econometrics*. Wiley.
#'
#' @seealso
#' [GCElm()], [summary.GCElm()]
#'
#' @examples
#' \donttest{
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
#' }
"_PACKAGE"
