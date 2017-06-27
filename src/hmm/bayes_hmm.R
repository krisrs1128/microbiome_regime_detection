#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Sample a Bayesian HMM using a block sampler. The advantage of this approach
## is that it lets us put a "sticky" prior on the transition probabilities.
##
## author: sankaran.kris@gmail.com
## date: 06/27/2017

## ---- libraries ----
library("mvtnorm")
library("MCMCpack")

## ---- sampling ----
#' @examples
#' K <- 20
#' gamma <- 2
#' alpha <- 0.9
#' kappa <- 3
#'
#' z <- c(rep(1, 10), rep(2, 5), rep(1, 10), rep(3, 3), rep(2, 10))
#' lambda <- list(zeta = .02, theta = c(0, 0), nu = 10, Delta = diag(c(0.1, 0.1)))
#' theta <- emission_parameters(K, lambda)
#' y <- emissions(z, theta)
#' plot(y[, 1], col = z)
#' z_clust <- kmeans(y, 20)$cluster
#' plot(y, col = 'white', asp = 1)
#' text(y, labels = z_clust)
#'
#' res <- block_sampler(y)
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

  n <- transition_counts(z, modes)
  Pi <- matrix(nrow = K, ncol = K, dimnames = list(modes, modes))

  for (k in modes) {
    alpha_new <- alpha + n[k, ]
    alpha_new[k] <- alpha_new[k] + kappa
    Pi[k, ] <- rdirichlet(1, alpha_new)[1, ]
  }

  Pi
}
