#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Plotting functions useful for hmm outputs.
##
## author: sankaran.kris@gmail.com
## date: 10/17/2017

gamma_mode <- function(gamma) {
  gamma_group <- gamma %>%
    group_by(ind, sample, asv) %>%
    summarise(k_max = K[which.max(gamma)])
  mu <- gamma %>%
    filter_(sprintf("sample == '%s'", gamma$sample[1])) %>%
    dplyr::select(K, mu) %>%
    unique()
  mu <- setNames(mu$mu, mu$K)

  gamma_group$mu <- mu[as.character(gamma_group$k_max)]
  gamma_group
}

melt_gamma <- function(gamma, dimn, samples, theta) {
  dimnames(gamma) <- dimn
  gamma <- gamma %>%
    melt(varnames = c("sample", "K", "asv"), value.name = "gamma") %>%
    as_data_frame() %>%
    left_join(samples)

  means <- data_frame(
    "K" = seq_along(theta),
    "mu" = sapply(theta, function(x) { x$mu })
  )
  k_order <- order(means$mu, decreasing = TRUE)

  gamma_mat <- gamma %>%
    dplyr::select(asv, sample, K, gamma) %>%
    unite(sample_K, sample, K) %>%
    spread(sample_K, gamma)

  hres <- hclust(dist(as.matrix(gamma_mat[, -1])))
  asv_order <- gamma_mat$asv[hres$order]

  gamma$mu <- means$mu[gamma$K]
  gamma$K <- factor(gamma$K, levels = k_order)
  gamma$asv <- factor(gamma$asv, levels = asv_order)
  gamma
}

extract_iteration_data <- function(samp_mcmc, n_iter, n_asv, n_sample) {
  K <- max(sapply(samp_mcmc, function(x) max(x$z)))
  z <- array(dim = c(n_sample, n_asv, n_iter))
  zgamma <- array(0, dim = c(n_sample, K, n_asv, n_iter))
  pi <- array(0, dim = c(K, K, n_iter))
  mu <- matrix(nrow = n_iter, ncol = K)
  sigma <- matrix(nrow = n_iter, ncol = K)
  for (iter in seq_len(n_iter)) {
    z[,, iter] <- samp_mcmc[[iter]]$z
    mu[iter, ] <- sapply(samp_mcmc[[iter]]$theta, function(x) { x$mu })
    sigma[iter, ] <- sapply(samp_mcmc[[iter]]$theta, function(x) { x$sigma })
    pi[,, iter] <- samp_mcmc[[iter]]$Pi
    for (k in seq_len(K)) {
      zgamma[, k,, iter] <- as.integer(z[,, iter] == k)
    }
  }

  list(
    "z" = z,
    "zgamma" = zgamma,
    "pi" = pi,
    "mu" = mu,
    "sigma" = sigma
  )
}

write_gif <- function(mz, name_base) {
  for (i in seq_len(max(mz$iter))) {
    cat(sprintf("saving iteration %s\n", i))
    p <- ggplot(mz %>% filter(iter == i)) +
      geom_tile(
        aes(x = sample, y = asv, fill = as.numeric(K))
      ) +
      scale_fill_viridis(option = "magma") +
      theme(
        axis.text = element_blank(),
        legend.position = "none"
      )
    ggsave(
      sprintf(
        "figure/%s_iter_%s.png",
        name_base,
        str_pad(i, 4, "left", "0")
      ),
      p, height = 9, width = 3.5
    )
  }
  system(
    sprintf(
      "convert -delay 15 -loop 0 figure/%s*.png figure/%s.gif",
      name_base, name_base
    )
  )
}

melt_z <- function(z, dimn, gamma) {
  colnames(z) <- dimn[[3]]
  rownames(z) <- dimn[[1]]
  mz <- melt(
    z,
    varnames = c("sample", "asv", "iter"),
    value.name = "K"
  ) %>%
    as_data_frame()
  mz$K <- factor(mz$K, levels(gamma$K))
  mz$asv <- factor(mz$asv, levels(gamma$asv))
  mz
}

extract_theta <- function(mu) {
  mu_df <- mu %>%
    melt(varnames = c("iter", "K"), value.name = "mu")
  theta <- mu_df %>%
    group_by(K) %>%
    summarise(mu = mean(mu))
  theta <- split(theta, seq(nrow(theta)))
  list(
    "mu_df" = mu_df,
    "theta" = theta
  )
}
