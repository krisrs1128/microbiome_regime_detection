#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Apply ordinary EM for HMMs to the antibiotics data set
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
res <- hmm_em(y, K, 20, lambda)
save(res, file = "hmm_em.rda")
