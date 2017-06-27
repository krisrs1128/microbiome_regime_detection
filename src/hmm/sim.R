#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Simulate toy data for experiments with different HMM models

simulate_data <- function() {
  z <- cbind(
    c(rep(1, 10), rep(2, 4), rep(1, 6), rep(3, 10)),
    c(rep(2, 6), rep(1, 5), rep(4, 5), rep(1, 4), rep(4, 10)),
    c(rep(4, 3), rep(1, 10), rep(2, 5), rep(4, 5), rep(3, 7))
  )[rep(1:30, each = 40), c(1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3)]
  time_len <- nrow(z)
  n <- ncol(z)
  K <- length(unique(z))

  ## parameters per state
  theta <- list(
    list("mu" = c(0, 0), "sigma" = diag(2)),
    list("mu" = c(-2, -1), "sigma" = diag(2)),
    list("mu" = c(2, 1), "sigma" = diag(2)),
    list("mu" = c(-1, -1), "sigma" = diag(2))
  )

  ## observed data
  y <- array(dim = c(time_len, n, 2))
  for (k in seq_along(theta)) {
    y[z == k] <- rmvnorm(
      sum(z == k),
      theta[[k]]$mu,
      theta[[k]]$sigma
    )
  }
  list("y" = y, "z" = z, "theta" = theta)
}
