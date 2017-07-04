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
source("em_hmm.R")
theme_set(ggscaffold::min_theme(list(
                        "legend_position" = "right",
                        "border_size" = 0.2
                      ))
          )
abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5 }, TRUE)

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
K <- 6
res <- hmm_em(y, K, 4, lambda)

###############################################################################
## inspect heatmap of states
###############################################################################
gamma <- res$gamma
dimnames(gamma) <- list(rownames(x), seq_len(K), colnames(x))
gamma <- gamma %>%
  melt(varnames = c("sample", "K", "rsv"), value.name = "gamma") %>%
  as_data_frame() %>%
  left_join(samples) %>%
  left_join(means)

means <- data_frame(
  "K" = 1:K,
  "mu" = sapply(res$theta, function(x) { x$mu })
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

cluster_cols <- c("#9cdea0", "#9cdeb6", "#9cdecc", "#9cdade", "#9cc4de", "#9caede")
ggplot(gamma) +
  geom_tile(
    aes(x = sample, y = rsv, alpha = gamma, fill = K)
  ) +
  scale_fill_manual(values = cluster_cols) +
  scale_alpha_continuous(range = c(0, 1)) +
  theme(axis.text = element_blank()) +
  facet_grid(K ~ ind, space = "free", scales = "free")

gamma_group <- gamma %>%
  group_by(ind, sample, rsv) %>%
  summarise(k_max = K[which.max(gamma)])
gamma_group$k_max <- factor(gamma_group$k_max, k_order)

ggplot(gamma_group) +
  geom_tile(
    aes(x = sample, y = rsv, fill = k_max)
  ) +
  scale_fill_manual(values = cluster_cols) +
  theme(axis.text = element_blank()) +
  facet_grid(. ~ ind, space = "free", scales = "free")

rownames(res$pi) <- 1:K
colnames(res$pi) <- 1:K
round(res$pi[k_order, k_order], 3)
