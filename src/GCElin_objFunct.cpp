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
NumericVector gceLinObj_eigen(SEXP lmb, SEXP env) {
        NumericVector lmbd(lmb);
        VectorXd lambda = as<Map<VectorXd> >(lmbd);
        Environment e = as<Environment>(env);
        MatrixXd X = as<MatrixXd>(e["X"]);
        VectorXd y = as<VectorXd>(e["y"]);
        MatrixXd Z = as<MatrixXd>(e["Z"]);
        MatrixXd v = as<MatrixXd>(e["v"]);
        double nu = as<double>(e["nu"]);
        MatrixXd p0 = as<MatrixXd>(e["p0"]);
        MatrixXd w0 = as<MatrixXd>(e["w0"]);
        int N = as<int>(e["N"]);
        int K = as<int>(e["K"]);
        // Omega
        MatrixXd lX = lambda.transpose()*X;
        int K1 = Z.rows();
        int s1 = Z.cols();
        MatrixXd numerator = MatrixXd(K1, s1).setZero();
        for (int k=0; k < K1; ++k) {
                double lXk = lX.col(k)[0];
                MatrixXd tmp1 = (Z.row(k).array()*lXk/(1-nu)).array().exp();
                numerator.row(k) = p0.row(k).array()*tmp1.array();
        } 
        MatrixXd Omega = numerator.rowwise().sum();
        MatrixXd Omega_sum = Omega.col(0).array().log().colwise().sum();
        // Psi
        MatrixXd lv = lambda*v.transpose();
        int N2 = lv.rows();
        int s2 = lv.cols();
        MatrixXd num_w = MatrixXd(N2,s2).setZero();
        for (int s=0; s < s2; ++s) {
                MatrixXd tmp2 = (lv.col(s)/nu).array().exp();
                num_w.col(s) = w0.col(s).array()*tmp2.array();
        } 
        MatrixXd Psi = num_w.rowwise().sum();
        MatrixXd Psi_sum = Psi.col(0).array().log().colwise().sum();
        // Likelihood
        MatrixXd lhs = -lambda.transpose()*y;
        MatrixXd MinValFun(lhs+(1-nu)*Omega_sum+nu*Psi_sum);
        NumericVector ret = as<NumericVector>(wrap(MinValFun));


//  Rcpp::Rcout<< "MinValFun=" << MinValFun << std::endl;
//  Rcpp::Rcout<< "tmp1=" << numerator.topLeftCorner(4,4) << std::endl;
//  Rcpp::Rcout<< "X=" << X(19999,1) << std::endl;
//  Rcpp::Rcout<< "Z=" << Z(2,1) << std::endl;
//  Rcpp::Rcout<< "v=" << v(0,0) << std::endl;
//  Rcpp::Rcout<< "K1=" << K1 << std::endl;
//  Rcpp::Rcout<< "s1=" << s1 << std::endl;
//  Rcpp::Rcout<< "lmbd=" << lmbd(0) << std::endl;
//  Rcpp::Rcout<< "lambda=" << lambda(0,0) << std::endl;
//  Rcpp::Rcout<< "lX=" << lX(0,0) << std::endl;

        return ret;
}


// definition
//' @export
// [[Rcpp::export(rng = false)]]
SEXP GCElin_objFunct() {
        BEGIN_RCPP
        typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
        return(XPtr<funcPtr>(new funcPtr(&gceLinObj_eigen)));
        END_RCPP
}
