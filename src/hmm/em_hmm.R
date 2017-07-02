#! /usr/bin/env Rscript

## File description ------------------------------------------------------------
## Fit a an HMM to parallel time series that share the same latent states, using
## the EM algorithm.
##
## author: sankaran.kris@gmail.com
## date: 07/02/2017

source("utils.R")
library("Rcpp")
sourceCpp("messages.cpp")

################################################################################
## Helpers that aren't key to the logic of the alborithm
################################################################################

normalize_log_space <- function(log_x) {
  log_x - lse(log_x)
}

normalize_rows_log <- function(log_x) {
  res <- matrix(0, nrow(log_x), ncol(log_x))
  for (i in seq_len(nrow(log_x))) {
    res[i, ] <- normalize_log_space(log_x[i, ])
  }
  res
}

merge_default_lambda <- function(lambda = list()) {
  default_lambda <- list(
    k0 = 0.1,
    m0 = rep(0, 2),
    nu0 = 3, ## p + 1 in our applications
    s0 = diag(1/(4 ^(1/2)) * 2, nrow = 2) # 1/k ^ (1/p) * column_sigma_hat
  )
  modifyList(default_lambda, lambda)
}

###############################################################################
## EM Algorithm
###############################################################################

#' EM Algorithm for Parallel HMM Sequences
#'
#' This is an implementation of standard EM for multiple sequencs each following
#' an HMM with a shared gaussian emission structure. We place prior on the
#' emission paramebers so the M-step is doing MAP estimation (see
#' merge_lambda_opts() above). The only big issue is that I'm not calculating
#' the loglikelihood at each iteration, which would be useful in comparing
#' solutions coming from different initializations, but it seems to work on the
#' toy examples I've created so far.
#'
#' @param y [matrix] A time x sequence numeric matrix
#' @param K [integer] The number of states
#' @param n_iter [integer] The number of iterations to run EM
#' @param lambda [list] A list of prior parameter options, see
#'   merge_lamba_opts()
#' @return res [list] A list with the following elements
#'   $gamma [array]: A time x state x sequence array giving the responsibility
#'   of cluster k for sequence j at time i
#'   $theta [list]: A list whose elements correspond to the estimated emission
#'     parameters for each state
#'   $pi [matrix]: The estimated transition probabilities between states
#'
#' @examples
#' sim <- simulate_data()
#' res <- hmm_em(sim$y)
#' plot(sim$y[, 1,], col = sim$z[, 1])
#' plot(sim$y[, 1 ,], col = apply(res$gamma[,, 1], 1, which.max))
#' plot(sim$y[, 1, 1], col = sim$z[,1])
#' plot(sim$y[, 1 , 1], col = apply(res$gamma[,, 1], 1, which.max))
#' image(sim$z)
#' image(apply(res$gamma, c(1, 3), which.max))
hmm_em <- function(y, K = 4, n_iter = 10, lambda = list()) {
  time_len <- nrow(y)
  n <- ncol(y)

  init <- initialize_states(y, K)
  theta <- init$theta
  pi <- init$n / rowSums(init$n)
  p0 <- setNames(rep(1 / K, K), seq_len(K))

  for (iter in seq_len(n_iter)) {
    cat(sprintf("iteration %s\n", iter))

    log_lik <- log_likelihood(y, theta)
    log_alpha <- array(dim = c(time_len, K, n))
    log_beta <- array(dim = c(time_len, K, n))
    log_gamma <- array(dim = c(time_len, K, n))
    log_xi <- array(dim = c(time_len - 1, K, K, n))

    ## E-step
    for (i in seq_len(n)) {
      log_alpha[,, i] <- forwards(pi, log_lik[,, i], p0)
      log_beta[,, i] <- backwards(pi, log_lik[,, i])
      log_gamma[,, i] <- normalize_rows_log(log_alpha[,, i] + log_beta[,, i])
      log_xi[,,, i] <- two_step_marginal(pi, log_lik[,, i], log_alpha[,, i], log_beta[,, i])
    }

    enjk <- expected_njk(log_xi)
    pi <- enjk / rowSums(enjk)
    p0 <- fitted_p0(as.matrix(log_gamma[1,, ]))
    theta <- expected_gaussian_param(y, exp(log_gamma), lambda)
  }

  list("theta" = theta, "gamma" = exp(log_gamma), "pi" = pi)
}

###############################################################################
## E-step: Infer latent states across all sequences in parallel
###############################################################################

#' @examples
#' pi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
#' lik <- matrix(runif(20), 10, 2)
#' forwards(pi, lik, c(0.5, 0.5))
forwards <- function(pi, log_lik, p0) {
  K <- nrow(pi)
  time_len <- nrow(log_lik)
  log_alpha <- matrix(0, time_len, K)

  ## base case
  log_alpha[1, ] <- log(p0) + log_lik[1, ]
  log_alpha[1, ] <- normalize_log_space(log_alpha[1, ])

  for (i in seq(2, time_len)) {
    log_alpha[i, ] <- log_lik[i, ] + log(t(pi) %*% exp(log_alpha[i - 1, ]))
    log_alpha[i, ] <- normalize_log_space(log_alpha[i, ])
  }

  log_alpha
}

#' @examples
#' pi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
#' lik <- matrix(runif(20), 10, 2)
#' backwards(pi, lik)
backwards <- function(pi, log_lik) {
  K <- nrow(pi)
  time_len <- nrow(log_lik)

  log_beta <- matrix(0, time_len, K)
  for (i in seq(time_len - 1, 1)) {
    for (k in seq_len(K)) {
      log_beta[i, k] <- lse(log_beta[i + 1, ] + log_lik[i + 1, ] + log(pi[k, ]))
    }
  }

  log_beta
}

#' @examples
#' pi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
#' lik <- matrix(runif(20), 10, 2)
#' beta <- exp(backwards(pi, lik))
#' alpha <- forwards(pi, lik)$alpha
#' two_step_marginal(pi, lik, alpha, beta)
#' @references Section 17.4.3.2 of "Machine Learning" by Murpy
two_step_marginal <- function(pi, log_lik, log_alpha, log_beta) {
  K <- nrow(pi)
  time_len <- nrow(log_lik)
  log_xi <- array(0, dim = c(time_len - 1, K, K))
  for (i in seq_len(time_len - 1)) {
    for (j in seq_len(K)) {
      for (k in seq_len(K)) {
        log_xi[i, j, k] <- log_alpha[i, k] + log_lik[i + 1, j] + log_beta[i + 1, j] + log(pi[k, j])
      }
    }
    log_xi[i,, ] <- normalize_log_space(log_xi[i,, ])
  }

  log_xi
}

fitted_p0 <- function(log_gamma1) {
  K <- nrow(log_gamma1)
  p0 <- vector(length = K)
  for (k in seq_len(K)) {
    p0[k] <- exp(lse(log_gamma1[k, ]))
  }
  p0
}

################################################################################
## M-step: Optimize emission parameters based on expected sufficient statistics.
################################################################################

log_likelihood <- function(y, theta) {
  K <- length(theta)
  time_len <- nrow(y)
  n <- ncol(y)

  log_lik <- array(dim = c(time_len, K, n))
  for (k in seq_len(K)) {
    for (i in seq_len(n)) {
      log_lik[, k, i] <- dmvn(y[, i,], theta[[k]]$mu, theta[[k]]$sigma, log = TRUE)
    }
  }
  log_lik
}

expected_njk <- function(log_xi) {
  K <- ncol(log_xi)
  log_njk <- matrix(nrow = K, ncol = K)
  for (j in seq_len(K)) {
    for (k in seq_len(K)) {
      log_njk[j, k] <- lse(log_xi[, j, k, ])
    }
  }

  exp(log_njk)
}

#' MAP estimation for gaussian parameters
#'
#' See section 11.4 in Murphy's Machine Learning textbook
expected_gaussian_param <- function(y, gamma, lambda = list()) {
  ## initialize prior parameters
  lambda <- merge_default_lambda(lambda)
  k0 <- lambda$k0
  m0 <- lambda$m0
  nu0 <- lambda$nu0
  s0 <- lambda$s0
  p <- dim(y)[3]
  K <- dim(gamma)[2]

  gamma <- aperm(gamma, c(1, 3, 2))
  y <- matrix(y, prod(dim(y)[1:2]), p)
  gamma <- matrix(gamma, prod(dim(gamma)[1:2]), K)

  ## each column is a weighted mean for cluster k
  cluster_weight <- colSums(gamma)
  y_means <- sweep(t(y) %*% gamma, 2, cluster_weight, "/")

  theta <- vector(mode = "list", length = K)
  for (k in seq_len(K)) {
    theta[[k]]$mu <- (cluster_weight[k] * y_means[, k] + k0 * m0) / (cluster_weight[k] + k0)

    sk <- 0
    for (i in seq_len(nrow(y))) {
      sk <- sk + gamma[i, k] * (y[i, ] - y_means[, k]) %*% t(y[i, ] - y_means[, k])
    }

    offset_coef <- (k0 * cluster_weight[k]) / (k0 + cluster_weight[k])
    theta[[k]]$sigma <- (1 / (nu0 + cluster_weight[k] + p + 2)) *
      (s0 + sk + offset_coef * (y_means[, k] - m0) %*% t(y_means[, k] - m0))
  }

  theta
}
