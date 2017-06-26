#! /usr/bin/env Rscript

## File description ------------------------------------------------------------
## Fit a an HMM to parallel time series that share the same latent states, using
## the EM algorithm.

################################################################################
## E-step: Compute smoothing probabilities and expected sufficient statistics
## via forwards backwards algorithm.
################################################################################

normalize <- function(x, return_sum = FALSE) {
  z <- sum(x)
  if (return_sum) {
    return (list("z" = z, "x" = x / z))
  }

  x / z
}

#' @examples
#' pi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
#' lik <- matrix(runif(20), 10, 2)
#' forwards(pi, lik)
forwards <- function(pi, lik, p0 = NULL) {
  K <- nrow(pi)
  time_len <- nrow(lik)
  if (is.null(p0)) {
    p0 <- rep(1 / K, K)
  }

  alpha <- matrix(0, time_len, K)
  Z <- vector(length = time_len)

  ## base case
  alpha[1, ] <- p0 * lik[1, ]
  Z[1] <- sum(alpha[1, ])
  alpha[1, ] <- alpha[1, ] / Z[1]

  for (i in seq(2, time_len)) {
    alpha[i, ] <- lik[i, ] * (t(pi) %*% alpha[i - 1, ])
    Z[i] <- sum(alpha[i, ])
    alpha[i, ] <- alpha[i, ] / Z[i]
  }

  list(alpha = alpha, Z = Z)
}

#' @examples
#' pi <- matrix(c(0.25, 0.75, 0.75, 0.25), 2)
#' lik <- matrix(runif(20), 10, 2)
#' backwards(pi, lik)
backwards <- function(pi, lik) {
  K <- nrow(pi)
  time_len <- nrow(lik)

  log_beta <- matrix(0, time_len, K)
  for (i in seq(time_len - 1, 1)) {
    log_beta[i, ] <- log(pi %*% (lik[i, ] * exp(log_beta[i + 1])))
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
two_step_marginal <- function(pi, lik, alpha, beta) {
  K <- nrow(pi)
  time_len <- nrow(lik)
  xi <- array(0, dim = c(time_len - 1, K, K))
  Z <- vector(length = time_len - 1)
  for (i in seq_len(time_len - 1)) {
    xi[i,, ] <- pi * (alpha[i, ] %*% t(lik[i + 1, ] * beta[i + 1, ]))
    Z[i] <- sum(xi[i,,])
    xi[i,, ] <- xi[i,,] / Z[i]
  }

  list("xi" = xi, "Z" = Z)
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

expected_njk <- function(xi) {
  apply(xi, c(2, 3), sum) # sum over times and samples
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
    theta[[k]]$sigma <- (yy_sums[k] - ns[k] * (theta[[k]]$mu ^ 2)) / ns[k]
  }

  theta
}

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

hmm_em <- function(Y, K = 5, n_iter = 10) {
  time_len <- nrow(Y)
  n <- ncol(Y)

  init <- initialize_states(Y, K)
  theta <- init$theta
  pi <- init$n / rowSums(init$n)
  p0 <- rep(1 / K, K)

  for (iter in seq_len(n_iter)) {
    log_lik <- log_likelihood(Y, theta)
    alpha <- array(dim = c(time_len, K, n))
    beta <- array(dim = c(time_len, K, n))
    gamma <- array(dim = c(time_len, K, n))
    xi <- array(dim = c(time_len - 1, K, K, n))

    ## E-step
    for (i in seq_len(n)) {
      alpha[,, i] <- forwards(pi, exp(log_lik[,, i]), p0)$alpha
      beta[,, i] <- exp(backwards(pi, exp(log_lik[,, i])))
      gamma[,, i] <- apply(alpha[,, i] * beta[,, i], 1, normalize)

      xi[,,, i] <- two_step_marginal(
        pi,
        exp(log_lik[,, i]),
        alpha[,, i],
        beta[,, i]
      )$xi
    }

    pi <- apply(expected_njk(xi), 1, normalize)
    p0 <- normalize(rowSums(gamma[1,, ]))
    theta <- expected_gaussian_param(Y, gamma)
  }

  list(
    "theta" = theta,
    "gamma" = gamma,
    "pi" = pi
  )
}

