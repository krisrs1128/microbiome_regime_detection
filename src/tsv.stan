## File description -------------------------------------------------------------
## Fit a simplified variant of the time series clustering model described in
## https://arxiv.org/abs/1505.01164

data {
  int<lower=1> K; // num topics
  int<lower=1> n; // num census tracts
  int<lower=1> T; // num timepoints
  vector[T] x[n]; // census tract dynamics
}

parameters {
  real lambda[n, K]; // latent factor loadings
  vector[T] eta[K]; // latent factor scores
  real sigma0;
  real mu_lambda;
  real sigma_lambda;
  simplex[K] theta[n]; // instead of z, we consider mixing proportions

  vector[n] a; // markov dynamics
}

model {
  for (k in 1:K) {
    eta[k] ~ normal(mu_lambda, sigma_lambda);
  }

  for (t in 1:(T - 1)) {
    for (i in 1:n) {
      for (k in 1:K) {
        target += theta[i, k] * normal_lpdf(x[i][t + 1] | a[i] * x[i][t] + lambda[i, k] * eta[k][t], sigma0);
      }
    }
  }
}
