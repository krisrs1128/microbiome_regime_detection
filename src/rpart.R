#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## This script is an experiment in using regression trees to create a
## simultaneous time x bacteria partitioning in the antibiotics data. The idea
## is to identify different dynamic regimes (normal vs. antibiotic) and taxa
## (bacteroides respond differently to firmicutes) without those labels being
## provided a priori.
##
## author: sankaran.kris@gmail.com
## date: 07/06/2017

###############################################################################
## Libraries and functions used below
###############################################################################
library("caret")
library("tidyverse")
library("viridis")
library("phyloseq")
library("dendextend")
library("vegan")
source("data.R")
theme_set(ggscaffold::min_theme())

pred_grid <- function(X, model, type = "raw") {
  x_grid <- expand.grid(
    "ind" = unique(X$ind),
    "leaf_ix" = unique(X$leaf_ix),
    "time" = unique(X$time)
  )
  y_hat <- predict(model, newdata = x_grid, type = type)
  if (type == "prob") {
    y_hat <- y_hat[, 2]
  }

  cbind(x_grid, y_hat)
}

plot_grid <- function(X, model, ...) {
  ggplot(pred_grid(X, model, ...)) +
    geom_tile(
      aes(x = time, y = leaf_ix, fill = y_hat)
    ) +
    scale_fill_viridis(option = "magma") +
    facet_grid(. ~ ind, scale = "free", space = "free") +
    theme(
      axis.text = element_blank()
    ) +
    scale_y_continuous(expand = c(0, 0))
}

###############################################################################
## Load and prepare the data of interest
###############################################################################
download.file(
  "https://github.com/krisrs1128/treelapse/raw/master/data/abt.rda",
  "../data/abt.rda"
)
abt <- get(load("../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5}, prune = TRUE)
x <- t(get_taxa(abt))
n <- nsamples(abt)

samples <- sample_data(abt) %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  as_data_frame()
mx <- melted_counts(x) %>%
  left_join(samples)

x_bin <- x > 0
D_euclidean <- dist(t(x), method = "euclidean")
D_jaccard <- vegdist(t(x_bin), method = "jaccard")
mix_tree <- hclust(0.5 * D_euclidean + 0.5 * D_jaccard)
mix_dendro <- reorder(as.dendrogram(mix_tree), -colMeans(x))

leaf_ix <- data_frame(
  "rsv" = mix_tree$labels,
  "leaf_ix" = order.dendrogram(mix_dendro)
)

mx  <- mx %>%
  left_join(leaf_ix)

###############################################################################
## train an rpart model across RSVs (mapped to their hclust index), subject, and
## time
###############################################################################
X <- mx %>%
  select(ind, time, leaf_ix)
y <- mx %>%
  .[["scaled"]]

train_opts <- list(
  "x" = X,
  "y" = y,
  "method" = "rpart",
  "tuneGrid" = data.frame("cp" = 1e-5),
  "trControl" = trainControl(method = "boot", number = 1)
)
rpart_model <- do.call(train, train_opts)
plot_grid(X, rpart_model$final_model)

###############################################################################
## different complexity penalization
###############################################################################
train_opts$tuneGrid <- data.frame("cp" = c(1e-4))
rpart_model <- do.call(train, train_opts)
plot_grid(X, rpart_model)

## notice differential recovery in F
train_opts$tuneGrid <- data.frame("cp" = c(5e-4))
rpart_model <- do.call(train, train_opts)
plot_grid(X, rpart_model)

train_opts$tuneGrid <- data.frame("cp" = c(7.5e-4))
rpart_model <- do.call(train, train_opts)
plot_grid(X, rpart_model)

###############################################################################
## Consider the conditionally positive models
###############################################################################
train_opts$tuneGrid <- data.frame("cp" = c(3e-4))
train_opts$x <- mx %>%
  filter(scaled > 0) %>%
  select(ind, time, leaf_ix)
train_opts$y <- mx %>%
  filter(scaled > 0) %>%
  .[["scaled"]]

## notice the lack of antibiotic timepoint splitting
rpart_model <- do.call(train, train_opts)
plot_grid(X, rpart_model)

###############################################################################
## Consider the binarized models
###############################################################################
train_opts$tuneGrid <- data.frame("cp" = c(3e-4))
train_opts$x <- mx %>%
  select(ind, time, leaf_ix)
train_opts$y <- mx %>%
  .[["present"]] %>%
  as.factor()
rpart_model <- do.call(train, train_opts)
plot_grid(X, rpart_model, type = "prob")

## notice that there are some that go from being not present to being present
## during the antibiotics time course
train_opts$tuneGrid <- data.frame("cp" = c(1e-4))
rpart_model <- do.call(train, train_opts)
plot_grid(X, rpart_model, type = "prob")
