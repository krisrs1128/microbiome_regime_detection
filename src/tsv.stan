## File description -------------------------------------------------------------
## Fit a simplified variant of the time series clustering model described in
## https://arxiv.org/abs/1505.01164

data {
  int<lower=1> K; // num topics
  int<lower=1> n; // num census tracts
  int<lower=1> T; // num timepoints
  vector[n] x[T]; // census tract dynamics
}

parameters {
  vector[n] lambda[K]; // latent factor loadings
  vector[T] eta[K]; // latent factor scores
  real sigma0;
  real mu_lambda;
  real sigma_lambda;
  matrix[n, K] beta; // instead of z, we consider mixing proportions

  vector[n] a; // markov dynamics
}

transformed parameters {
  matrix[n, K] theta;
  for (i in 1:n) {
    theta[i] = inv_logit(beta[i]);
  }
}

model {
  for (k in 1:K) {
    eta[k] ~ normal(mu_lambda, sigma_lambda);
  }

  for (t in 1:(T - 1)) {
    for (k in 1:K) {
      target += theta[, k] * normal_lpdf(x[t + 1] | a .* x[t] + eta[k][t] * lambda[k], sigma0);
    }
  }
}
