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
library("viridis")
library("jsonlite")
library("stringr")
theme_set(ggscaffold::min_theme(list(
                        "legend_position" = "right",
                        "border_size" = 0.2
                      ))
          )

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
        aes(x = sample, y = rsv, fill = as.numeric(K))
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
samples <- sample_data(abt) %>%
  data.frame() %>%
  rownames_to_column("sample")

res <- get(load("hmm_em.rda"))
gamma <- res$gamma
K <- nrow(res$pi)
dimn <- list(sample_names(abt), seq_len(K), taxa_names(abt))
gamma <- melt_gamma(gamma, dimn, samples, res$theta)

p <- ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = as.numeric(K))
  ) +
  scale_fill_viridis(option = "magma") +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank(), legend.position = "none") +
  facet_grid(K ~ ind, space = "free", scales = "free")
ggsave("~/Desktop/lab_meetings/20170705/figure/hmm_probs.pdf", height = 3, width = 2)

p <- ggplot(gamma_mode(gamma)) +
  geom_tile(
    aes(x = sample, y = rsv, fill = as.numeric(k_max))
  ) +
  scale_fill_viridis(option = "magma") +
  theme(axis.text = element_blank(), legend.position = "none") +
  facet_grid(. ~ ind, space = "free", scales = "free")
ggsave("~/Desktop/lab_meetings/20170705/figure/hmm_mode.pdf", height = 3, width = 2)

rownames(res$pi) <- 1:K
colnames(res$pi) <- 1:K
round(res$pi[levels(gamma$K), levels(gamma$K)], 3)

###############################################################################
## Block sampler for bayesian HMM
###############################################################################
samp_mcmc <- readLines("bayes_kappa_4.txt")
samp_mcmc <- samp_mcmc[seq(1, length(samp_mcmc), 1)] %>%
  lapply(fromJSON)
K <- nrow(samp_mcmc[[1]]$Pi)
dimn <- list(sample_names(abt), seq_len(K), taxa_names(abt))

samp_data <- extract_iteration_data(
  samp_mcmc,
  length(samp_mcmc),
  ntaxa(abt),
  nsamples(abt)
)
rm(samp_mcmc)
gamma <- apply(samp_data$zgamma, c(1, 2, 3), mean)
samp_data$zgamma <- NULL

theta <- extract_theta(samp_data$mu)
gamma <- melt_gamma(gamma, dimn, samples, theta$theta)

theta$mu_df$K <- factor(theta$mu_df$K, levels(gamma$K))
ggplot(theta$mu_df) +
  geom_histogram(
    aes(x = mu, fill = K),
    binwidth = 0.01
  ) +
  facet_wrap(~ K, scales = "free")

pi <- samp_data$pi %>%
  melt(
    varnames = c("i", "j", "iter"),
    value.name = "pi_ij"
  )
pi$i <- factor(pi$i, levels(gamma$K))
pi$j <- factor(pi$j, levels(gamma$K))

pi_df <- pi %>%
  group_by(i, j) %>%
  summarise(pi_mean = mean(pi_ij), se = sd(pi_ij))

pi_mean <- pi_df %>%
  select(-se) %>%
  spread(j, pi_mean) %>%
  ungroup %>%
  select(-i)
pi_se <- pi_df %>%
  select(-pi_mean) %>%
  spread(j, se) %>%
  ungroup %>%
  select(-i)

round(pi_mean, 3)
round(pi_se, 3)

ggplot(pi) +
  geom_histogram(aes(x = pi_ij), binwidth = 0.005) +
  facet_grid(i ~ j)

ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = K)
  ) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(K ~ ind, space = "free", scales = "free")

p <- ggplot(gamma_mode(gamma)) +
  geom_tile(
    aes(x = sample, y = rsv, fill = as.numeric(k_max))
  ) +
  scale_fill_viridis(option = "magma") +
  theme(axis.text = element_blank(), legend.position = "none") +
  facet_grid(. ~ ind, space = "free", scales = "free")
ggsave("~/Desktop/lab_meetings/20170705/figure/bayes_mode.pdf", height = 3, width = 2)

mz <- melt_z(samp_data$z[,, 1000:1050], dimn, gamma)
p <- ggplot(mz %>% filter(rsv %in% levels(gamma$rsv)[sample(1:600, 4)])) +
  geom_tile(
    aes(x = sample, y = iter, fill = as.numeric(K))
  ) +
  scale_fill_viridis(option = "magma") +
  facet_wrap(~rsv) +
  theme(axis.text = element_blank())
ggsave("~/Desktop/lab_meetings/20170705/figure/gibbs_samples.pdf", height = 1.8, width = 3)
write_gif(mz, "bayes_hmm_z")

###############################################################################
## Figures for HDP-HMM results
###############################################################################
samp_mcmc <- readLines("hdp_kappa_01.txt") %>%
  lapply(fromJSON)
K <- nrow(samp_mcmc[[1]]$Pi)
dimn <- list(sample_names(abt), seq_len(K), taxa_names(abt))

samp_data <- extract_iteration_data(
  samp_mcmc,
  length(samp_mcmc),
  ntaxa(abt),
  nsamples(abt)
)
rm(samp_mcmc)
gamma <- apply(samp_data$zgamma, c(1, 2, 3), mean)
samp_data$zgamma <- NULL

theta <- extract_theta(samp_data$mu)
dimn[[2]] <- as.character(seq_len(ncol(gamma)))
gamma <- melt_gamma(gamma, dimn, samples, theta$theta)

theta$mu_df$K <- factor(theta$mu_df$K, levels(gamma$K))
ggplot(theta$mu_df) +
  geom_histogram(
    aes(x = mu, fill = K),
  ) +
  facet_wrap(~ K)

ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = K)
  ) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(K ~ ind, space = "free", scales = "free")

p <- ggplot(gamma_mode(gamma)) +
  geom_tile(
    aes(x = sample, y = rsv, fill = as.numeric(k_max))
  ) +
  scale_fill_viridis(option = "magma") +
  theme(axis.text = element_blank(), legend.position = "none") +
  facet_grid(. ~ ind, space = "free", scales = "free")
ggsave("~/Desktop/lab_meetings/20170705/figure/hdp_antibiotics_mode.pdf", height = 3, width = 2)

## study mixing in z
mz <- melt_z(samp_data$z[,, 1000:1050], dimn, gamma)
ggplot(mz %>% filter(rsv %in% levels(gamma$rsv)[1:3])) +
  geom_tile(
    aes(x = sample, y = iter, fill = as.numeric(K))
  ) +
  scale_fill_viridis(option = "magma") +
  facet_wrap(~rsv) +
  theme(axis.text = element_blank(),
        legend.position = "none")
write_gif(mz, "hdp_hmm_z")

## look at transition probabilities
pi <- samp_data$pi %>%
  melt(
    varnames = c("i", "j", "iter"),
    value.name = "pi_ij"
  )
pi$i <- factor(pi$i, levels(gamma$K))
pi$j <- factor(pi$j, levels(gamma$K))

pi_df <- pi %>%
  group_by(i, j) %>%
  summarise(pi_mean = mean(pi_ij), se = sd(pi_ij))

pi_mean <- pi_df %>%
  select(-se) %>%
  spread(j, pi_mean) %>%
  ungroup %>%
  select(-i)
pi_se <- pi_df %>%
  select(-pi_mean) %>%
  spread(j, se) %>%
  ungroup %>%
  select(-i)

round(pi_mean, 3)
round(pi_se, 3)
