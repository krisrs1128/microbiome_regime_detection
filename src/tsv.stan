## File description -------------------------------------------------------------
## Fit a simplified variant of the time series clustering model described in
## https://arxiv.org/abs/1505.01164

data {
  int<lower=1> K; // num topics
  int<lower=1> n; // num census tracts
  int<lower=1> T; // num timepoints
  real X[T, n]; // census tract dynamics
}

parameters {
  real eps[T, n]; // clustered innovations
  real lambda[n, K]; // latent factor loadings
  real[T] eta[K]; // latent factor scores
  real sigma0;
  real mu_lambda;
  real sigma_lambda;
}
