#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Apply the different HMM algorithms to the antibiotics data set.
##
## author: sankaran.kris@gmail.com
## date: 10/18/2017

###############################################################################
## setup
###############################################################################
library("tidyverse")
library("reshape2")
library("phyloseq")
library("viridis")
library("jsonlite")
library("forcats")
library("stringr")
source("plot.R")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0.04, unit = "cm"),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

###############################################################################
## inspect heatmap of states
###############################################################################
abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5 }, TRUE)
samples <- sample_data(abt) %>%
  data.frame() %>%
  mutate(time = str_pad(time, 2, "left", 0)) %>%
  unite(sample, ind, time, sep = "", remove = FALSE)

res <- get(load("hmm_em.rda"))
gamma <- res$gamma
K <- nrow(res$pi)
dimn <- list(samples$sample, seq_len(K), taxa_names(abt))
gamma <- melt_gamma(gamma, dimn, samples, res$theta)

viri_scale <- scale_fill_viridis(
  option = "magma",
  direction = -1,
  limits = c(0, 5.2),
  breaks = c(0, 2, 4)
)
legend_guide <- guides(
  fill = guide_colorbar(
    barheight = 0.2,
    ticks = FALSE
  ),
  alpha = guide_legend(
    keywidth = unit(0.3, "cm"),
    ticks = FALSE
  )
)

ggplot(gamma) +
  geom_tile(
    aes(x = asv, y = sample, alpha = gamma, fill = mu)
  ) +
  viri_scale +
  legend_guide +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(
    axis.text = element_blank(),
    legend.position = "bottom",
    panel.border = element_rect(size = 1, fill = "transparent")
  ) +
  facet_grid(ind ~ K, space = "free", scales = "free")
ggsave(
  "../../doc//figure/hmm_probs.png",
  width = 6.67,
  height = 3.01
)

taxa <- tax_table(abt) %>%
  data.frame() %>%
  rownames_to_column("asv")
taxa$asv <- factor(taxa$asv, levels(gamma$asv))
taxa$family <- taxa$Taxon_5
taxa$family[taxa$family == ""] <- NA
taxa$family <- fct_lump(taxa$family, 7, ties.method = "first")
taxa$family <- taxa$family %>%
  recode(
    Alcaligenaceae_Sutterella = "Sutterella",
    Peptostreptococcaceae_1 = "Peptostreptococcaceae"
  )
taxa$family[is.na(taxa$family)] <- "Other"
taxa$family <- factor(
  taxa$family,
  names(sort(table(taxa$family), decreasing = TRUE))
)

gm <- gamma_mode(gamma) %>%
  left_join(taxa)
ggplot(gm) +
  geom_tile(
    aes(x = asv, y = sample, fill = mu)
  ) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      col = family
    ),
    fill = "transparent",
    size = 0.9
  ) +
  viri_scale +
  legend_guide +
  theme(
    axis.text = element_blank(),
    legend.position = "bottom",
    panel.spacing.x = unit(0.1 ,"cm"),
    panel.border = element_blank(),
    strip.text.x = element_blank()
  ) +
  facet_grid(ind ~ family, space = "free", scales = "free")

ggsave(
  "../../doc/figure/hmm_mode.png",
  width = 6.87,
  height = 3.53
)

k_order <- extract_theta(res$theta)$mu_df %>%
  na.omit() %>%
  select(L1, mu) %>%
  arrange(mu) %>%
  .[["L1"]]

res$pi <- res$pi[k_order, k_order]
rownames(res$pi) <- 1:K
colnames(res$pi) <- 1:K
round(res$pi, 3)

###############################################################################
## Block sampler for bayesian HMM
###############################################################################
file_pipe <- pipe("awk 'BEGIN{i=0}{i++;if (i%10==2) print $1}' < bayes_kappa_4.txt")
samp_mcmc <- readLines(file_pipe)
samp_mcmc <- lapply(samp_mcmc, fromJSON)
K <- nrow(samp_mcmc[[1]]$Pi)
dimn <- list(samples$sample, seq_len(K), taxa_names(abt))

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
  dplyr::select(-se) %>%
  spread(j, pi_mean) %>%
  ungroup
pi_se <- pi_df %>%
  select(-pi_mean) %>%
  spread(j, se) %>%
  ungroup

round(pi_mean[, -1], 3)
round(pi_se[, -1], 3)

ggplot(pi) +
  geom_histogram(aes(x = pi_ij), binwidth = 0.005) +
  facet_grid(i ~ j)

ggplot(gamma) +
  geom_tile(
    aes(x = asv, y = sample, alpha = gamma, fill = K)
  ) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(ind ~ K, space = "free", scales = "free")

levels(taxa$asv) <- levels(gamma$asv)
gm <- gamma_mode(gamma) %>%
  left_join(taxa)
ggplot(gm) +
  geom_tile(
    aes(x = asv, y = sample, fill = mu)
  ) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      col = family
    ),
    fill = "transparent",
    size = 0.9
  ) +
  viri_scale +
  legend_guide +
  theme(
    axis.text = element_blank(),
    legend.position = "bottom",
    panel.spacing.x = unit(0.1 ,"cm"),
    strip.text.x = element_blank()
  ) +
  facet_grid(ind ~ family, space = "free", scales = "free")
ggsave(
  "../../doc/figure/bayes_mode.png",
  width = 6.87,
  height = 3.53
)

mz <- melt_z(samp_data$z, dimn, gamma)
ggplot(mz %>% filter(asv %in% levels(gamma$asv)[sample(1:600, 9)])) +
  geom_tile(
    aes(x = sample, y = iter, fill = K)
  ) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~asv) +
  theme(axis.text = element_blank())
ggsave(
  "../../doc//figure/bayes_gibbs_samples.png",
  width = 8.26,
  height = 4.151
)
## write_gif(mz, "bayes_hmm_z")

###############################################################################
## Figures for HDP-HMM results
###############################################################################
file_pipe <- pipe("awk 'BEGIN{i=0}{i++;if (i%10==2) print $1}' < hdp_kappa_01.txt")
samp_mcmc <- readLines(file_pipe) %>%
  lapply(fromJSON)
K <- nrow(samp_mcmc[[1]]$Pi)
dimn <- list(samples$sample, seq_len(K), taxa_names(abt))

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
  ) +
  facet_wrap(~ K)

ggplot(gamma) +
  geom_tile(
    aes(x = asv, y = sample, alpha = gamma, fill = K)
  ) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(ind ~ K, space = "free", scales = "free")

levels(taxa$asv) <- levels(gamma$asv)
gm <- gamma_mode(gamma) %>%
  left_join(taxa)
ggplot(gm) +
  geom_tile(
    aes(x = asv, y = sample, fill = mu)
  ) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      col = family
    ),
    fill = "transparent",
    size = 0.9
  ) +
  viri_scale +
  legend_guide +
  theme(
    axis.text = element_blank(),
    legend.position = "bottom",
    panel.spacing.x = unit(0.1 ,"cm"),
    strip.text.x = element_blank()
  ) +
  facet_grid(ind ~ family, space = "free", scales = "free")

ggsave(
  "../../doc/figure/hdp_mode.png",
  width = 6.87,
  height = 3.53
)

## study mixing in z
mz <- melt_z(samp_data$z, dimn, gamma)
ggplot(mz %>% filter(asv %in% levels(gamma$asv)[sample(1:600, 9)])) +
  geom_tile(
    aes(x = sample, y = iter, fill = K)
  ) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~asv) +
  theme(axis.text = element_blank())
ggsave(
  "../../doc//figure/hdp_gibbs_samples.png",
  width = 8.26,
  height = 4.151
)
## write_gif(mz, "hdp_hmm_z")

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
  dplyr::summarise(pi_mean = mean(pi_ij), se = sd(pi_ij))

pi_mean <- pi_df %>%
  dplyr::select(-se) %>%
  spread(j, pi_mean)
pi_se <- pi_df %>%
  dplyr::select(-pi_mean) %>%
  spread(j, se)

round(pi_mean[, -1], 3)
round(pi_se[, -1], 3)
