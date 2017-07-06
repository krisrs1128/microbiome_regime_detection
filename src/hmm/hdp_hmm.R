#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Block sampler for the HDP-HMM. Based on Algorithm 10 of ``Bayesian
## Nonparametric Learning of Complex Dynamical Phenomena''
##
## author: sankaran.kris@gmail.com
## date: 6/21/2017

## ---- libraries ----
library("mvnfast")
library("jsonlite")
library("MCMCpack")
library("Rcpp")
sourceCpp("messages.cpp")
source("utils.R")

################################################################################
## Default options for hyperparameters and prior
################################################################################

merge_default_hyper <- function(opts = list()) {
  default_opts <- list(
    "L" = 20,
    "n_iter" = 20,
    "theta_iter" = 2,
    "alpha" = 1,
    "gamma" = 2,
    "kappa" = 3,
    "outpath" = file.path(getwd(), "samples.txt")
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

###############################################################################
## Non-auxiliary variable sampling functions
###############################################################################

#'  Block sampling for the Sticky HDP-HMM
#'
#' This is an implementation of Emily Fox's block sampling algorithm for the
#' Sticky HDP-HMM. Note that inference is taking place across many sequences
#' simultaneously, which is a minor change from the original paper -- all that
#' changes is that the z-sampling step is performed across each separately.
#'
#' @examples
#' sim <- simulate_data()
#' hyper <- list(
#'   n_iter = 200,
#'   L = 10,
#'   gamma = 1e-4, ## controls total number of states
#'   alpha = 1e-2, ## controls number of states any one state can transition to
#'   kappa = 0.01
#' )
#' res <- block_sampler(sim$y[, 1:20, ], hyper)
#' plot(sim$y[,1, 1], col = res$z[, 1])
#' plot(sim$y[,1, ], col = res$z[, 1], asp = 1)
#' image(res$z)
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
  beta <- setNames(rep(0.1, L), 1:L)

  for (iter in seq_len(hyper$n_iter)) {
    cat(sprintf("iteration %s/%s\n", iter, hyper$n_iter))

    for (i in seq_len(ncol(y))) {
      msg <- messages(Pi, as.matrix(y[, i,, drop = FALSE]), theta)
      z[, i] <- sample_z(Pi, as.matrix(y[, i,, drop = FALSE]), theta, msg)
    }

    m <- sample_m(z, hyper$alpha, beta, hyper$kappa)
    w <- sample_override(
      diag(m),
      hyper$kappa / (hyper$kappa + hyper$alpha),
      beta
    )
    beta <- sample_beta(m, w, hyper$gamma)

    Pi <- sample_pi(z, hyper$alpha, beta, hyper$kappa)
    theta <- sample_theta(y, z, theta, lambda, hyper$theta_iter)
    state <- list("z" = z, "beta" = beta, "theta" = theta, "Pi" = Pi)

    cat(
      sprintf("%s\n", toJSON(c("iter" = iter, state))),
      file = hyper$outpath,
      append = iter != 1
    )
  }

  state
}

sample_beta <- function(m, w, gamma) {
  m_bar <- m
  diag(m_bar) <- diag(m_bar) - w

  setNames(
    rdirichlet(1, gamma / ncol(m_bar) + colSums(m_bar))[1, ],
    colnames(m)
  )
}

sample_pi <- function(z, alpha, beta, kappa) {
  modes <- names(beta)
  n <- 0
  for (i in seq_len(ncol(z))) {
    n <- n + transition_counts(z[, i], modes)
  }

  Pi <- matrix(0, length(modes), length(modes),
               dimnames = list(modes, modes))
  for (l in modes) {
    u <- alpha * beta + n[l, ]
    u[l] <- u[l] + kappa
    Pi[l, ] <- rdirichlet(1, u)[1, ]
  }
  Pi
}

###############################################################################
## Functions used to sample the auxiliary variables in the HDP-HMM. These are in
## their own file because this step is shared in both the direct assignment and
## the block samplers.
###############################################################################

#' Random CRT distributed variable
#'
#' This is the distribution of the number of tables associated with a CRP after
#' m people have been seated.
#'
#' See https://dukespace.lib.duke.edu/dspace/handle/10161/7204
rcrt <- function(gamma, m) {
  probs <- gamma / (seq_len(m) - 1 + gamma)
  sum(rbinom(m, 1, probs))
}

sample_m <- function(z, alpha, beta, kappa) {
  modes <- names(beta)
  m <- matrix(0, length(modes), length(modes),
              dimnames = list(modes, modes))
  n <- 0
  for (i in seq_len(ncol(z))) {
    n <- n + transition_counts(z[, i], modes)
  }

  for (j in modes) {
    for (k in modes) {
      m[j, k] <- rcrt(alpha * beta[k] + kappa * (j == k), n[j, k])
    }
  }

  m
}

sample_override <- function(m_diag, rho, beta) {
  K <- length(m_diag)
  w <- vector(length = K)
  for (k in seq_len(K)) {
    w[k] <- rbinom(1, m_diag[k], rho / (rho + beta[k] * (1 - rho)))
  }

  w
}
