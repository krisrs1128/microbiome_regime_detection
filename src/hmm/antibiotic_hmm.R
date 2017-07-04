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
abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5 }, TRUE)
K <- 6
cluster_cols <- c("#9cdea0", "#9cdeb6", "#9cdecc", "#9cdade", "#9cc4de", "#9caede")

###############################################################################
## Some utilities
###############################################################################

melted_counts <- function(x) {
  x %>% data.frame() %>%
    rownames_to_column("sample") %>%
    as_data_frame %>%
    gather(key = "rsv", value = "value", -sample) %>%
    mutate(
      scaled = asinh(value),
      present = ifelse(value > 0, 1, 0)
    )
}

taxa_labels <- function(taxa) {
  taxa <- taxa %>%
    data.frame() %>%
    rownames_to_column("rsv")
  taxa$label <- taxa$Taxon_5
  taxa$label[is.na(taxa$Taxon_5)] <- taxa$Taxon_4[is.na(taxa$Taxon_5)]
  ordered_labels <- names(sort(table(taxa$label), decreasing = TRUE))
  taxa_levels <- c(ordered_labels, "other")
  taxa$label <- factor(taxa$label, taxa_levels)
  taxa$label[!c(taxa$label %in% ordered_labels[1:6])] <- "other"
  taxa$label[taxa$label == ""] <- "other"
  taxa
}

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
    "K" = 1:K,
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

write_gif <- function(mz) {
  for (i in seq_len(max(mz$iter))) {
    cat(sprintf("saving iteration %s\n", i))
    p <- ggplot(mz %>% filter(iter == i)) +
      geom_tile(
        aes(x = sample, y = rsv, fill = K)
      ) +
      scale_fill_manual(values = cluster_cols) +
      theme(axis.text = element_blank())
    ggsave(
      sprintf("figure/z_iter_%s.png", str_pad(i, 4, "left", "0")),
      p, height = 9, width = 4
    )
  }
  system("convert -delay 10 -loop 0 figure/*.png figure/z.gif")
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
## Extract data and fit HMM to ordinary data
###############################################################################
x_df <- abt %>%
  get_taxa() %>%
  t() %>%
  melted_counts()

samples <- sample_data(abt) %>%
  data.frame() %>%
  rownames_to_column("sample")

x_df <- x_df %>%
  select(sample, rsv, scaled) %>%
  spread(rsv, scaled) %>%
  left_join(samples) %>%
  arrange(ind, time)

sample_names <- x_df$sample
x <- x_df %>%
  select(-sample, -ind, -time, -condition) %>%
  as.matrix()
rownames(x) <- sample_names

y <- array(x, dim = c(dim(x), 1))
lambda <- list("mu" = mean(x), "nu0" = 2, "s0" = 1, "m0" = 0)
#res <- hmm_em(y, K, 4, lambda)

###############################################################################
## inspect heatmap of states
###############################################################################
gamma <- res$gamma
dimn <- list(rownames(x), seq_len(K), colnames(x))
gamma <- melt_gamma(gamma, dimn, samples, res$theta)

ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = K)
  ) +
  scale_fill_manual(values = cluster_cols) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(K ~ ind, space = "free", scales = "free")

gamma_group <- gamma_mode(gamma)

ggplot(gamma_group) +
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
samp_mcmc <- readLines("bayes_kappa_4.txt") %>%
  lapply(fromJSON)

samp_data <- extract_iteration_data(
  samp_mcmc,
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

gamma_group <- gamma_mode(gamma)

ggplot(gamma_group) +
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
write_gif(mz)
