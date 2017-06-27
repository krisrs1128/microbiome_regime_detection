#! /usr/bin/env Rscript

## File description ------------------------------------------------------------
## Fit a an HMM to parallel time series that share the same latent states, using
## the EM algorithm.

################################################################################
## E-step: Compute smoothing probabilities and expected sufficient statistics
## via forwards backwards algorithm.
################################################################################

## ---- initialization ----
#' @examples
#' z <- c(1, 1, 2, 1, 1, 3, 3, 3, 1)
#' transition_counts(z)
transition_counts <- function(z) {
  modes <- sort(unique(z))
  n <- matrix(0, nrow = length(modes), ncol = length(modes),
              dimnames = list(modes, modes))
  time_len <- length(z)

  z <- as.character(z)
  for (i in seq_len(time_len - 1)) {
    n[z[i], z[i + 1]] <- n[z[i], z[i + 1]] + 1
  }
  n
}

initialize_states <- function(Y, K) {
  y <- c(Y)
  y_clust <- kmeans(y, K)

  theta <- vector(mode = "list", length = K)
  for (k in seq_len(K)) {
    theta[[k]]$mu <- mean(y[y_clust$cluster == k])
    theta[[k]]$sigma <- sd(y[y_clust$cluster == k])
  }

  list(
    "theta" = theta,
    "n" = transition_counts(y_clust$cluster)
  )
}

## ---- e-step ----
lse <- function(log_x) {
  m <- max(log_x)
  m + log(sum(exp(log_x - m)))
}

normalize_log_space <- function(log_x) {
  log_x - lse(log_x)
}

normalize <- function(x) {
  x / sum(x)
}

normalize_rows <- function(x) {
  x / rowSums(x)
}

normalize_rows_log <- function(log_x) {
  res <- matrix(0, nrow(log_x), ncol(log_x))
  for (i in seq_len(nrow(log_x))) {
    res[i, ] <- normalize_log_space(log_x[i, ])
  }
  res
}

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

log_likelihood <- function(Y, theta) {
  K <- length(theta)
  time_len <- nrow(Y)
  n <- ncol(Y)

  log_lik <- array(dim = c(time_len, K, n))
  for (k in seq_len(K)) {
    log_lik[, k, ] <- dnorm(Y, theta[[k]]$mu, theta[[k]]$sigma, log = TRUE)
  }
  log_lik
}

expected_njk <- function(log_xi) {
  K <- ncol(log_xi)
  log_njk <- matrix(nrow = K, ncol = K)
  for (j in seq_len(K)) {
    for (k in seq_len(K)) {
      log_njk[j, k] <- lse(log_xi[, j, k,])
    }
  }

  exp(log_njk)
}

expected_nj <- function(alpha, beta) {
  gamma <- alpha * beta
  apply(gamma, 2, sum) # sum over times and samples
}

expected_gaussian_param <- function(Y, gamma) {
  K <- ncol(gamma)
  ns <- vector(length = K)
  y_sums <- vector(length = K)
  yy_sums <- vector(length = K)
  for (k in seq_len(K)) {
    y_sums[k] <- sum(gamma[, k, ] * Y)
    yy_sums[k] <- sum(gamma[, k, ] * (Y ^ 2))
    ns[k] <- sum(gamma[, k, ])
}

  theta <- vector(mode = "list", length = k)
  for (k in seq_len(K)) {
    theta[[k]]$mu <- y_sums[k] / ns[k]
    theta[[k]]$sigma <- sqrt((yy_sums[k] - ns[k] * (theta[[k]]$mu ^ 2)) / ns[k])
  }

  theta
}

hmm_em <- function(Y, K = 4, n_iter = 10) {
  time_len <- nrow(Y)
  n <- ncol(Y)

  init <- initialize_states(Y, K)
  theta <- init$theta
  pi <- (init$n) / rowSums(init$n)
  p0 <- setNames(rep(1 / K, K), seq_len(K))

  for (iter in seq_len(n_iter)) {
    cat(sprintf("iteration %s\n", iter))

    log_lik <- log_likelihood(Y, theta)
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

    pi <- normalize_rows(expected_njk(log_xi))
    p0 <- fitted_p0(as.matrix(log_gamma[1,, ]))
    theta <- expected_gaussian_param(Y, log_gamma)
  }

  list("theta" = theta, "gamma" = exp(log_gamma), "pi" = pi)
}
