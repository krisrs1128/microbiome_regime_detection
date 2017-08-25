#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Apply the different HMM algorithms to the antibiotics data set.
##
## author: sankaran.kris@gmail.com
## date: 07/04/2017

###############################################################################
## setup
###############################################################################
source("bayes_hmm.R")
library("reshape2")
library("phyloseq")
library("jsonlite")
library("tidyverse")
set.seed(705)
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5 }, TRUE)
K <- 4

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

###############################################################################
## Extract data
###############################################################################
x_df <- abt %>%
  get_taxa() %>%
  t() %>%
  melted_counts()

samples <- sample_data(abt) %>%
  data.frame() %>%
  rownames_to_column("sample")

x_df <- x_df %>%
  dplyr::select(sample, rsv, scaled) %>%
  spread(rsv, scaled) %>%
  left_join(samples) %>%
  arrange(ind, time)

sample_names <- x_df$sample
x <- x_df %>%
  dplyr::select(-sample, -ind, -time, -condition) %>%
  as.matrix()
rownames(x) <- sample_names

y <- array(x, dim = c(dim(x), 1))

###############################################################################
## Block sampler for bayesian HMM
###############################################################################
hyper <- list(
  "L" = K,
  "n_iter" = 2000,
  "alpha" = setNames(rep(1, K), seq_len(K)),
  "kappa" = 4,
  "outpath" = "bayes_kappa_4.txt"
)
lambda <- list("mu0" = mean(y), "sigma0" = sd(y), "nu" = 2, "delta" = matrix(1))
res <- block_sampler(y, hyper, lambda)
