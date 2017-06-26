#! /usr/bin/env Rscript

## File description ------------------------------------------------------------
## Fit a an HMM to parallel time series that share the same latent states, using
## the EM algorithm.

################################################################################
## E-step: Compute smoothing probabilities and expected sufficient statistics
## via forwards backwards algorithm.
################################################################################

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
    alpha[i, ] <- lik[i, ] * (pi %*% alpha[i - 1, ])
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

################################################################################
## M-step: Optimize emission parameters based on expected sufficient statistics.
################################################################################
