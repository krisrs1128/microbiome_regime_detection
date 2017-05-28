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
  vector[T] eps[n]; // clustered innovations
  real lambda[n, K]; // latent factor loadings
  vector[T] eta[K]; // latent factor scores
  real sigma0;
  real mu_lambda;
  real sigma_lambda;
  simplex[K] theta[n]; // instead of z, we consider mixing proportions
}

model {
  for (i in 1:n) {
    for (k in 1:K) {
      target += theta[i, k] * normal_lpdf(eps[i] | lambda[i, k] * eta[k], sigma0)
    }
  }
}
