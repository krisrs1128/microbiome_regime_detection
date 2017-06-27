// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
#include <vector>
#include <random>

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double lse(vec log_x) {
  double m = max(log_x);
  return m + log(sum(exp(log_x - m)));
}

//' http://gallery.rcpp.org/articles/dmvnorm_arma/
// [[Rcpp::export]]
double dmvnrm(vec x, vec mean, mat sigma) {
  int xdim = x.n_elem;
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  vec z = rooti * (x - mean) ;
  double out = constants - 0.5 * sum(z%z) + rootisum;
  return(out);
}

// [[Rcpp::export]]
vec multi_dmvnorm(vec yt, Rcpp::List theta) {
  int L = theta.length();
  vec y_dens = zeros<vec>(L);
  for (int l = 0; l < L; l++) {
    Rcpp::List theta_l = theta[l];
    vec mu = theta_l["mu"];
    mat sigma = theta_l["sigma"];
    y_dens[l] = dmvnrm(yt, mu, sigma);
  }

  return y_dens;
}

vec row_(mat x, int i) {
  return conv_to<vec>::from(x.row(i));
}

// [[Rcpp::export]]
mat messages(mat Pi, mat y, Rcpp::List theta) {
  int time_len = y.n_rows;
  int K = Pi.n_cols;
  mat log_msg = zeros<mat>(time_len, K);
  mat logPi = log(Pi);

  for (int i = time_len - 2; i > 0; i--) {
    vec log_y_dens = multi_dmvnorm(row_(y, i), theta);
    for (int k = 0; k < K; k++) {
      log_msg(i, k) = lse(row_(logPi, k) + log_y_dens + row_(log_msg, i + 1));
    }
  }

  return log_msg;
}

//[[Rcpp::export]]
vec sample_z(mat Pi, mat y, Rcpp::List theta, mat msg) {
  int time_len = y.n_rows;
  int K = Pi.n_cols;
  vec z = ones<vec>(time_len);
  mat logPi = log(Pi);

  std::random_device rd;
  std::mt19937 gen(rd());


  for (int i = 2; i < time_len; i++) {
    vec log_y_dens = multi_dmvnorm(row_(y, i), theta);
    vec log_f = row_(logPi, z(i - 1)) + log_y_dens + row_(msg, i);
    std::vector<double> p = conv_to<std::vector<double> >::from(exp(log_f - lse(log_f)));
    std::discrete_distribution<> d(p.begin(), p.end());

    z[i] = d(gen);
  }

  // convert since R is 1 indexed
  for (int i = 0; i < time_len; i++) {
    z[i] += 1;
  }

  return z;
}
