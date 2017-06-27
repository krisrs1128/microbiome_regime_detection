
## ---- initialization ----
#' @examples
#' z <- c(1, 1, 2, 1, 1, 3, 3, 3, 1)
#' transition_counts(z)
transition_counts <- function(z, modes = NULL) {
  if (is.null(modes)) {
    modes <- sort(unique(z))
  }

  n <- matrix(0, nrow = length(modes), ncol = length(modes),
              dimnames = list(modes, modes))
  time_len <- length(z)

  z <- as.character(z)
  for (i in seq_len(time_len - 1)) {
    n[z[i], z[i + 1]] <- n[z[i], z[i + 1]] + 1
  }
  n
}

initialize_states <- function(y, K) {
  y_clust <- kmeans(y, K)
  theta <- vector(mode = "list", length = K)
  for (k in seq_len(K)) {
    theta[[k]]$mu <- mean(y[y_clust$cluster == k, ])
    theta[[k]]$sigma <- 0.75 * cov(y[y_clust$cluster == k, ])
  }

  list(
    "z": y_clust$cluster,
    "theta" = theta,
    "n" = transition_counts(y_clust$cluster)
  )
}

## ---- sampling-latent-states ----
messages <- function(Pi, y, theta) {
  time_len <- nrow(y)
  modes <- colnames(Pi)
  log_msg <- matrix(0, time_len, nrow(Pi),
                    dimnames = list(1:time_len, modes))

  for (i in seq(time_len - 1, 1)) {
    log_y_dens <- multi_dmvnorm(y[i, ], theta)
    log_msg[i, ] <- log(Pi %*% exp(log_y_dens + log_msg[i + 1, ]))
  }
  log_msg
}

sample_z <- function(Pi, y, theta, msg) {
  time_len <- nrow(y)
  z <- rep(1, time_len)
  for (i in seq(2, time_len)) {
    log_y_dens <- multi_dmvnorm(y[i, ], theta)
    f <- exp(log(Pi[z[i - 1], ]) + log_y_dens + msg[i, ])
    z[i] <- sample(seq_along(f), 1, prob = f / sum(f))
  }
  z
}

## ---- sampling-posterior-parameters ----
sample_mu <- function(mu0, sigma0, y, sigma) {
  if (nrow(y) == 0) {
    return (rmvnorm(1, mu0, sigma0)[1, ])
  }

  sigma_bar <- solve(solve(sigma0) + nrow(y) * solve(sigma))
  mu_bar <- sigma_bar %*% (solve(sigma0) %*% mu0 + solve(sigma) %*% colSums(y))
  rmvnorm(1, mu_bar, sigma_bar)[1, ]
}

sample_sigma <- function(nu, delta, y, mu) {
  if (nrow(y) == 0) {
    return (riwish(nu, nu * delta))
  }

  nu_bar <- nu + nrow(y)
  e <- y - rep(1, nrow(y)) %*% matrix(mu, nrow = 1)
  nu_delta_bar <- nu * delta + t(e) %*% e
  riwish(nu_bar, nu_delta_bar)
}

sample_theta <- function(y, z, theta, lambda, n_iter) {
  modes <- names(theta)
  z <- as.character(z)

  for (i in seq_len(n_iter)) {
    for (l in modes) {
      theta[[l]]$mu <- sample_mu(
        lambda$mu0,
        lambda$sigma0,
        y[z == l,, drop = FALSE],
        theta[[l]]$sigma
      )

      theta[[l]]$sigma <- sample_sigma(
        lambda$nu,
        lambda$delta,
        y[z == l,, drop = FALSE],
        theta[[l]]$mu
      )
    }
  }

  theta
}

## ---- probability-densities ----
multi_dmvnorm <- function(yt, theta) {
  modes <- names(theta)
  y_dens <- setNames(seq_along(modes), modes)
  for (l in modes) {
    y_dens[l] <- mvtnorm::dmtvnorm(
                            yt,
                            theta[[l]]$mu,
                            theta[[l]]$sigma,
                            log = TRUE
                          )
  }
  y_dens
}
