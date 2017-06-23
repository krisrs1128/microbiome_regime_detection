#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## This script is an experiment in using regression trees to create a
## simultaneous time x bacteria partitioning in the antibiotics data. The idea
## is to identify different dynamic regimes (normal vs. antibiotic) and taxa
## (bacteroides respond differently to firmicutes) without those labels being
## provided a priori.
##
## author: sankaran.kris@gmail.com
## date: 06/22/2017

## ---- libraries ----
library("caret")
library("tidyverse")
library("phyloseq")
source("data.R")

## ---- data ----
abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 70}, prune = TRUE) %>%
  transform_sample_counts(asinh)
x <- t(get_taxa(abt))
n <- nsamples(abt)

samples <- sample_data(abt) %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  as_data_frame()
mx <- melted_counts(x) %>%
  left_join(samples)

phylo_mds <- cmdscale(1 - cophenetic(phy_tree(abt))) %>%
  as.data.frame() %>%
  rownames_to_column("rsv") %>%
  as_data_frame() %>%
  rename(phylo_ix = V1)

mx  <- mx %>%
  left_join(phylo_mds)

## ---- model ----
X <- mx %>%
  filter(value > 0) %>%
  select(time, rsv)
y <- mx %>%
  filter(value > 0) %>%
  .[["value"]]

tune_grid <- data.frame("cp" = c(0.01))
rpart_model <- train(X, y, method = "rpart", tuneGrid = tune_grid)
plot(rpart_model$finalModel)

X <- mx %>%
  filter(value > 0) %>%
  select(time, phylo_ix)
rpart_model <- train(X, y, method = "rpart", tuneGrid = tune_grid)
