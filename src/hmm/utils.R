
## ---- initialization ----
#' @examples
#' z <- c(1, 1, 2, 1, 1, 3, 3, 3, 1)
#' transition_counts(z)
transition_counts <- function(z, modes = NULL) {
  if (is.null(modes)) {
    modes <- sort(unique(z))
  }

  z <- factor(z, levels = modes)
  table(head(z, -1), tail(z, -1))
}

initialize_states <- function(y, L) {
  y0 <- y
  y <- matrix(y, prod(dim(y)[1:2]), dim(y)[3]) # unfold array to matrix

  y_clust <- kmeans(y, L)
  theta <- setNames(vector(mode = "list", length = L), 1:L)
  for (l in seq_len(L)) {
    theta[[l]]$mu <- colMeans(y[y_clust$cluster == l,, drop = FALSE])
    theta[[l]]$sigma <- 0.75 * cov(y[y_clust$cluster == l,, drop = FALSE])
  }

  list(
    "z" = matrix(y_clust$cluster, nrow(y0), ncol(y0)),
    "theta" = theta,
    "n" = transition_counts(y_clust$cluster)
  )
}

## ---- sampling-posterior-parameters ----
sample_mu <- function(mu0, sigma0, y, sigma) {
  if (nrow(y) == 0) {
    return (rmvn(1, mu0, sigma0)[1, ])
  }

  sigma_bar <- solve(solve(sigma0) + nrow(y) * solve(sigma))
  mu_bar <- sigma_bar %*% (solve(sigma0) %*% mu0 + solve(sigma) %*% colSums(y))
  rmvn(1, mu_bar, sigma_bar)[1, ]
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
  y <- matrix(y, prod(dim(y)[1:2]), dim(y)[3])
  z <- as.character(z) # unwraps into character vector

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

## ---- toy-example ----
simulate_data <- function() {
  z <- cbind(
    c(rep(1, 10), rep(2, 4), rep(1, 6), rep(3, 10)),
    c(rep(2, 6), rep(1, 5), rep(4, 5), rep(1, 4), rep(4, 10)),
    c(rep(4, 3), rep(1, 10), rep(2, 5), rep(4, 5), rep(3, 7))
  )[rep(1:30, each = 4), c(rep(1, 20), rep(2, 10), rep(3, 40))]
  time_len <- nrow(z)
  n <- ncol(z)
  K <- length(unique(z))

  ## parameters per state
  theta <- list(
    list("mu" = c(-2, -1), "sigma" = diag(2)),
    list("mu" = c(-1, -1), "sigma" = diag(2)),
    list("mu" = c(0, 0), "sigma" = diag(2)),
    list("mu" = c(2, 1), "sigma" = diag(2))
  )

  ## observed data
  y <- array(dim = c(time_len, n, 2))
  for (k in seq_along(theta)) {
    y[z == k] <- rmvn(
      sum(z == k),
      theta[[k]]$mu,
      theta[[k]]$sigma
    )
  }
  list("y" = y, "z" = z, "theta" = theta)
}
