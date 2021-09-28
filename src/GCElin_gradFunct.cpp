// includes from the plugin
// #include <Rcpp.h>
  // [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;
  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;


// user includes
NumericVector gceLinGrad_eigen(SEXP lmb, SEXP env) {

  NumericVector lmbd(lmb);
  VectorXd lambda = as<Map<VectorXd> >(lmbd);
  Environment ev = as<Environment>(env);
  MatrixXd X = as<MatrixXd>(ev["X"]);
  VectorXd y = as<VectorXd>(ev["y"]);
  MatrixXd Z = as<MatrixXd>(ev["Z"]);
  MatrixXd v = as<MatrixXd>(ev["v"]);
  double nu = as<double>(ev["nu"]);
  MatrixXd p0 = as<MatrixXd>(ev["p0"]);
  MatrixXd w0 = as<MatrixXd>(ev["w0"]);
  int N = as<int>(ev["N"]);
  int K = as<int>(ev["K"]);
  int M = as<int>(ev["M"]);
  int J = as<int>(ev["J"]);
  // Omega
  MatrixXd lX = lambda.transpose()*X;
  int K1 = Z.rows();
  int m1 = Z.cols();
  MatrixXd numerator = MatrixXd(K1, m1).setZero();
  for (int k=0; k < K1; ++k) {
    double lXk = lX.col(k)[0];
    MatrixXd tmp1 = (Z.row(k).array()*lXk/(1-nu)).array().exp();
    numerator.row(k) = p0.row(k).array()*tmp1.array();
  } 
  MatrixXd Omega = numerator.rowwise().sum();
  // Matrice p
  MatrixXd p = MatrixXd(K1, m1).setZero();
  for (int m=0; m < m1; ++m) {
    MatrixXd tmp2 = numerator.col(m).array()/Omega.col(0).array();
    p.col(m) = tmp2;
  }
  // beta
  MatrixXd beta = (Z.array()*p.array()).rowwise().sum();
  // Psi
  MatrixXd lv = lambda*v.transpose();
  int N2 = lv.rows();
  int m2 = lv.cols();
  MatrixXd num_w = MatrixXd(N2,m2).setZero();
  for (int m=0; m < m2; ++m) {
    MatrixXd tmp3 = (lv.col(m)/nu).array().exp();
    num_w.col(m) = w0.col(m).array()*tmp3.array();
  } 
  MatrixXd Psi = num_w.rowwise().sum();
  // Matrice w ed e
  MatrixXd w = MatrixXd(N2,m2).setZero();
  for (int m=0; m < m2; ++m) {
    MatrixXd tmp4 = num_w.col(m).array()/Psi.col(0).array();
    w.col(m) = tmp4;
  }
  MatrixXd e =  (w*v).rowwise().sum();
  VectorXd gradMinValFun = -y + X*beta + e;
  NumericVector ret = as<NumericVector>(wrap(gradMinValFun));
  return ret;
}

// definition
//' @export
// [[Rcpp::export(rng = false)]]
SEXP GCElin_gradFunct() {
  typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&gceLinGrad_eigen)));
}
