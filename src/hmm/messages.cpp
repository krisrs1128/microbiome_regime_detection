#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double lse(vec log_x) {
  double m = max(log_x);
  return m + log(sum(exp(log_x - m)));
}

// [[Rcpp::export]]
mat messages(mat Pi, mat y, Rcpp::List theta) {
  int time_len = y.n_rows;
  int K = Pi.n_cols;
  mat log_msg = zeros<mat>(time_len, K);

  for (int i = time_len - 1; i > 1; i--) {
    // vec log_y_dens = multi_dmvnorm(y.row(i), theta);
    vec log_y_dens = zeros<vec>(10);
      for (int k = 0; k < K; k++) {
        log_msg[i, k] = lse(log(Pi.row(k)) + log_y_dens + log_msg.row(i + 1));
      }
  }

  return log_msg;
}

