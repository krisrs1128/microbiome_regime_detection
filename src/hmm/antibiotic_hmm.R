#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Apply the different HMM algorithms to the antibiotics data set.
##
## author: sankaran.kris@gmail.com
## date: 07/04/2017

###############################################################################
## setup
###############################################################################
library("tidyverse")
library("reshape2")
library("phyloseq")
library("jsonlite")
library("stringr")
source("em_hmm.R")
theme_set(ggscaffold::min_theme(list(
                        "legend_position" = "right",
                        "border_size" = 0.2
                      ))
          )
col_fun <- colorRampPalette(c("#9cdea0", "#9caede"))

###############################################################################
## Some utilities
###############################################################################
gamma_mode <- function(gamma) {
  gamma_group <- gamma %>%
    group_by(ind, sample, rsv) %>%
    summarise(k_max = K[which.max(gamma)])
  gamma_group$k_max <- factor(gamma_group$k_max, levels(gamma$K))
  gamma_group
}

melt_gamma <- function(gamma, dimn, samples, theta) {
  dimnames(gamma) <- dimn
  gamma <- gamma %>%
    melt(varnames = c("sample", "K", "rsv"), value.name = "gamma") %>%
    as_data_frame() %>%
    left_join(samples)

  means <- data_frame(
    "K" = seq_along(theta),
    "mu" = sapply(theta, function(x) { x$mu })
  )
  k_order <- order(means$mu, decreasing = TRUE)

  gamma_mat <- gamma %>%
    select(rsv, sample, K, gamma) %>%
    unite(sample_K, sample, K) %>%
    spread(sample_K, gamma)

  hres <- hclust(dist(as.matrix(gamma_mat[, -1])))
  rsv_order <- gamma_mat$rsv[hres$order]

  gamma$K <- factor(gamma$K, levels = k_order)
  gamma$rsv <- factor(gamma$rsv, levels = rsv_order)
  gamma
}

extract_iteration_data <- function(samp_mcmc, n_iter, n_rsv, n_sample) {
  K <- max(sapply(samp_mcmc, function(x) max(x$z)))
  z <- array(dim = c(n_sample, n_rsv, n_iter))
  zgamma <- array(0, dim = c(n_sample, K, n_rsv, n_iter))
  mu <- matrix(nrow = n_iter, ncol = K)
  sigma <- matrix(nrow = n_iter, ncol = K)
  for (iter in seq_len(n_iter)) {
    z[,, iter] <- samp_mcmc[[iter]]$z
    mu[iter, ] <- sapply(samp_mcmc[[iter]]$theta, function(x) { x$mu })
    sigma[iter, ] <- sapply(samp_mcmc[[iter]]$theta, function(x) { x$sigma })

    for (k in seq_len(K)) {
      zgamma[, k,, iter] <- as.integer(z[,, iter] == k)
    }
  }

  list(
    "z" = z,
    "zgamma" = zgamma,
    "mu" = mu,
    "sigma" = sigma
  )
}

write_gif <- function(mz, name_base) {
  for (i in seq_len(max(mz$iter))) {
    cat(sprintf("saving iteration %s\n", i))
    p <- ggplot(mz %>% filter(iter == i)) +
      geom_tile(
        aes(x = sample, y = rsv, fill = K)
      ) +
      scale_fill_manual(values = cluster_cols) +
      theme(axis.text = element_blank())
    ggsave(
      sprintf(
        "figure/%s_iter_%s.png",
        name_base,
        str_pad(i, 4, "left", "0")
      ),
      p, height = 9, width = 4
    )
  }
  system(
    sprintf(
      "convert -delay 15 -loop 0 figure/%s*.png figure/%s.gif",
      name_base, name_base
    )
  )
}

melt_z <- function(z, x, gamma) {
  colnames(z) <- colnames(x)
  rownames(z) <- rownames(x)
  mz <- melt(
    z,
    varnames = c("sample", "rsv", "iter"),
    value.name = "K"
  ) %>%
    as_data_frame()
  mz$K <- factor(mz$K, levels(gamma$K))
  mz$rsv <- factor(mz$rsv, levels(gamma$rsv))
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

###############################################################################
## inspect heatmap of states
###############################################################################
abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5 }, TRUE)
res <- get(load("hmm_em.rda"))
gamma <- res$gamma
gamma <- melt_gamma(gamma, dimn, samples, res$theta)
K <- nrow(res$Pi)
dimn <- list(sample_names(abt), seq_len(K), taxa_names(abt))
cluster_cols <- col_fun(K)

ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = K)
  ) +
  scale_fill_manual(values = cluster_cols) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(K ~ ind, space = "free", scales = "free")

ggplot(gamma_mode(gamma)) +
  geom_tile(
    aes(x = sample, y = rsv, fill = k_max)
  ) +
  scale_fill_manual(values = cluster_cols) +
  theme(axis.text = element_blank()) +
  facet_grid(. ~ ind, space = "free", scales = "free")

rownames(res$pi) <- 1:K
colnames(res$pi) <- 1:K
round(res$pi[levels(gamma$K), levels(gamma$K)], 3)

###############################################################################
## Block sampler for bayesian HMM
###############################################################################
samp_mcmc <- readLines("bayes_kappa_4.txt")
samp_mcmc <- samp_mcmc[seq(1, length(samp_mcmc), 10)] %>%
  lapply(fromJSON)
K <- nrow(samp_mcmc[[1]]$Pi)
dimn <- list(sample_names(abt), seq_len(K), taxa_names(abt))
cluster_cols <- col_fun(K)

samp_data <- extract_iteration_data(
 rsamp_mcmc,
  length(samp_mcmc),
  ncol(y),
  nrow(y)
)
gamma <- apply(samp_data$zgamma, c(1, 2, 3), mean)

theta <- extract_theta(samp_data$mu)
gamma <- melt_gamma(gamma, dimn, samples, theta$theta)

theta$mu_df$K <- factor(theta$mu_df$K, levels(gamma$K))
ggplot(theta$mu_df) +
  geom_histogram(
    aes(x = mu, fill = K),
  ) +
  scale_fill_manual(values = cluster_cols) +
  facet_wrap(~ K)

ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = K)
  ) +
  scale_fill_manual(values = cluster_cols) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(K ~ ind, space = "free", scales = "free")

ggplot(gamma_mode(gamma)) +
  geom_tile(
    aes(x = sample, y = rsv, fill = k_max)
  ) +
  scale_fill_manual(values = cluster_cols) +
  theme(axis.text = element_blank()) +
  facet_grid(. ~ ind, space = "free", scales = "free")

mz <- melt_z(samp_data$z, x, gamma)
ggplot(mz %>% filter(rsv %in% levels(gamma$rsv)[1:3])) +
  geom_tile(
    aes(x = sample, y = iter, fill = K)
  ) +
  scale_fill_manual(values = cluster_cols) +
  facet_wrap(~rsv) +
  theme(axis.text = element_blank())
write_gif(mz, "bayes_hmm_z")

###############################################################################
## Figures for HDP-HMM results
###############################################################################
samp_mcmc <- readLines("hdp_kappa_01.txt") %>%
  lapply(fromJSON)
K <- nrow(samp_mcmc[[1]]$Pi)
dimn <- list(sample_names(abt), seq_len(K), taxa_names(abt))
cluster_cols <- col_fun(K)

samp_data <- extract_iteration_data(
  samp_mcmc,
  length(samp_mcmc),
  ncol(y),
  nrow(y)
)
gamma <- apply(samp_data$zgamma[,,, 100:200], c(1, 2, 3), mean)

theta <- extract_theta(samp_data$mu)
dimn[[2]] <- as.character(seq_len(ncol(gamma)))
gamma <- melt_gamma(gamma, dimn, samples, theta$theta)

theta$mu_df$K <- factor(theta$mu_df$K, levels(gamma$K))
cluster_cols <- col_fun(length(unique(gamma$K)))
ggplot(theta$mu_df) +
  geom_histogram(
    aes(x = mu, fill = K),
  ) +
  scale_fill_manual(values = cluster_cols) +
  facet_wrap(~ K)

ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = K)
  ) +
  scale_fill_manual(values = cluster_cols) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(K ~ ind, space = "free", scales = "free")

ggplot(gamma_mode(gamma)) +
  geom_tile(
    aes(x = sample, y = rsv, fill = k_max)
  ) +
  scale_fill_manual(values = cluster_cols) +
  theme(axis.text = element_blank()) +
  facet_grid(. ~ ind, space = "free", scales = "free")

mz <- melt_z(samp_data$z, x, gamma)
ggplot(mz %>% filter(rsv %in% levels(gamma$rsv)[1:3])) +
  geom_tile(
    aes(x = sample, y = iter, fill = K)
  ) +
  scale_fill_manual(values = cluster_cols) +
  facet_wrap(~rsv) +
  theme(axis.text = element_blank())
write_gif(mz, "hdp_hmm_z")
