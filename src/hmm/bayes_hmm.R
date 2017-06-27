#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Sample a Bayesian HMM using a block sampler. The advantage of this approach
## is that it lets us put a "sticky" prior on the transition probabilities.
##
## author: sankaran.kris@gmail.com
## date: 06/27/2017

## ---- libraries ----
library("mvnfast")
library("MCMCpack")
source("utils.R")
Rcpp::sourceCpp("messages.cpp")

## ---- sampling ----
#' @examples
#' sim <- simulate_data()
#' block_sampler(sim$y, list(K = 4, n_iter = 200))
block_sampler <- function(y, hyper = list(), lambda = list()) {
  ## merge default opts
  hyper <- merge_default_hyper(hyper)
  lambda <- merge_default_lambda(lambda)

  ## initialize state space
  L <- hyper$L
  init <- initialize_states(y, L)

  Pi <- init$n / rowSums(init$n)
  z <- init$z
  theta <- init$theta

  for (iter in seq_len(hyper$n_iter)) {
    cat(sprintf("iteration %s\n", iter))

    for (i in seq_len(ncol(y))) {
      msg <- messages(Pi, y[, i,], theta)
      z[, i] <- sample_z(Pi, y[, i,], theta, msg)
    }

    Pi <- sample_pi(z, hyper$alpha, hyper$kappa)
    theta <- sample_theta(y, z, theta, lambda, hyper$theta_iter)
    state <- list(z = z, beta = beta, theta = theta)
  }

  state
}

## ---- defaults ----
merge_default_hyper <- function(opts = list()) {
  default_opts <- list(
    "L" = 4,
    "n_iter" = 500,
    "theta_iter" = 2,
    "alpha" = setNames(rep(1, 4), seq_len(4)),
    "kappa" = 1
  )
  modifyList(default_opts, opts)
}

merge_default_lambda <- function(opts = list()) {
  default_opts <- list(
    "mu0" = c(0, 0),
    "sigma0" = diag(2),
    "nu" = 3,
    "delta" = matrix(c(1, 0, 0, 1), nrow = 2)
  )
  modifyList(default_opts, opts)
}

sample_pi <- function(z, alpha, kappa) {
  K <- length(alpha)
  modes <- names(alpha)

  n <- transition_counts(z[, 1], modes)
  for (i in seq(2, ncol(z))) {
    n <- n + transition_counts(z[, i], modes)
  }

  Pi <- matrix(nrow = K, ncol = K, dimnames = list(modes, modes))
  for (k in modes) {
    alpha_new <- alpha + n[k, ]
    alpha_new[k] <- alpha_new[k] + kappa
    Pi[k, ] <- rdirichlet(1, alpha_new)[1, ]
  }

  Pi
}
