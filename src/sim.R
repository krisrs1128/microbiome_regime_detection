#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulate data similar to that described in https://arxiv.org/abs/1505.01164
## Main difference is we don't have situations like "homes within a tract". So,
## we just model clustered census tracts.

## ---- libraries ----
library("plyr")
library("dplyr")
library("ggplot2")
library("tibble")
library("tidyr")
library("rstan")
theme_set(ggscaffold::min_theme())

## ---- utils ----
sim_x <- function(a, eps, x0 = NULL) {
  n <- nrow(eps)
  p <- ncol(eps)

  X <- matrix(0, n, p)
  if (is.null(x0)) {
    X[1, ] <- rnorm(p)
  }

  for (i in seq_len(n - 1)) {
    X[i + 1, ] <- a * X[i, ] + eps[i + 1, ]
  }

  X
}

sim_eps <- function(n, p, k, sigma0, mu_lambda, sigma_lambda) {
  lambda_sim <- sim_lambda_z(p, k, mu_lambda, sigma_lambda)

  eps <- matrix(0, n, p)
  for (i in seq_len(n)) {
    eps[i, ] <- lambda_sim$lambda_z %*% rnorm(k)
  }

  list(
    "eps" = eps + matrix(rnorm(n * p, 0, sigma0), n, p), ## add \tilde{eps} contribution
    "lambda" = lambda_sim$lambda,
    "z" = lambda_sim$z
  )
}

sim_lambda_z <- function(p, k, mu_lambda, sigma_lambda) {
  Lambda <- matrix(rnorm(p * k, mu_lambda, sigma_lambda), p, k)

  Z <- matrix(0, p, k)
  z <- sample(seq_len(k), p, replace = TRUE)
  for (i in seq_len(p)) {
    Z[i, z[i]] <- 1
  }

  list(
    "lambda_z" = Lambda * Z,
    "lambda" = Lambda,
    "z" = z
  )
}

## ---- parameters ----
## 100 times, 30 tracts, 4 tract-clusters
params <- list(
  n = 100,
  p = 200,
  k = 3,
  mu_lambda = 1,
  sigma_lambda = 0.8,
  sigma0 = 1,
  a = rbeta(200, 4, 1)
)

## ---- simulate-epsilon ----
eps_res <- do.call(sim_eps, params[setdiff(names(params), "a")])
mapping <- setNames(eps_res$z, paste0("V", seq_len(params$p)))
eps_df <- as_data_frame(eps_res$eps) %>%
  mutate(time = row_number()) %>%
  gather(key = tract, value = value, -time) %>%
  mutate(k = mapping[tract])

p <- ggplot(eps_df) +
  geom_line(
    aes(x = time, y = value, group = tract),
    size = 0.3, alpha = 0.7
  ) +
  facet_wrap(~ k)

image(cor(eps_res$eps[, order(eps_res$z)]))

## ---- simulate-x ----
X <- sim_x(params$a, eps_res$eps)

X_df <- as_data_frame(X) %>%
  mutate(time = row_number()) %>%
  gather(key = tract, value = value, -time) %>%
  mutate(k = mapping[tract])

p <- ggplot(X_df) +
  geom_line(
    aes(x = time, y = value, group = tract),
    size = 0.3, alpha = 0.7
  ) +
  facet_wrap(~ k)

## ---- fit-model ----
stan_data <- list(
  "K" = params$k,
  "n" = params$p,
  "T" = params$n,
  "x" = t(X),
  "alpha" = c(0.4, 0.4, 0.4)
)

fit <- vb(stan_model("ts_cluster.stan"), stan_data, iter = 1000, adapt_engaged = FALSE, eta = 1)

estimates <- extract(fit)
plot(params$a, colMeans(estimates$a))
theta_hat <- apply(estimates$theta, c(2, 3), mean)

table(
  truth = eps_res$z,
  estimates = apply(theta_hat, 1, which.max)
)

names(estimates)
mean(estimates$mu_lambda)
mean(estimates$sigma_lambda)
mean(estimates$sigma0)
