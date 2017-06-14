#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some experiments using wavelet denoising on the antibiotics data, in an
## attempt to cluster the timepoints across different bacteria.

## ---- libraries ----
library("phyloseq")
library("wavethresh")
library("reshape2") library("tidyverse")
library("forcats")
theme_set(ggscaffold::min_theme())

## ---- utils ----
interpolate_dyadic <- function(y) {
  n <- length(y)
  n_interp <- 2 ^ ceiling(log(n, 2))
  f_y <- approxfun(seq_along(y), y, yleft = 0, yright = 0)
  f_y(seq(1, n, length.out = n_interp))
}

stacked_detail <- function(ywd, min_lev = 1, max_lev = NULL, scale_lev = FALSE) {
  if (is.null(max_lev)) {
    max_lev <- ywd$nlevel
  }

  n_interp <- length(ywd$D) + 1
  stacked <- matrix(0, n_interp, max_lev - min_lev + 1)
  for (l in min_lev:(max_lev - 1)) {
    stacked[, l] <- rep(accessD(ywd, l), each = n_interp / 2 ^ l)
  }
  if (scale_lev) {
    stacked <- scale(stacked)
    stacked[is.na(stacked)] <- 0
  }

  stacked
}

detail_array <- function(y, ...) {
  n <- nrow(y)
  p <- ncol(y)
  D <- array(0, dim = c(2 ^ ceiling(log(n, 2)), p, ceiling(log(n, 2))))
  transforms <- list()
  for (j in seq_len(p)) {
    transforms[[j]] <- x[, j] %>%
      interpolate_dyadic %>%
      wd(filter.number = 1, family = "DaubExPhase")
    D[, j,] <- stacked_detail(transforms[[j]], ...)
  }
  list("D" = D, "transforms" = transforms)
}

cluster_detail_array <- function(D, K = 8) {
  mdetail <- D %>%
    melt(varnames = c("interp_time", "rsv", "coef_level"))

  cluster_data <- mdetail %>%
    mutate(coef_level = paste0("L", coef_level)) %>%
    spread(coef_level, value)

  cluster_mat <- cluster_data %>%
    select(starts_with("L")) %>%
    as.matrix

  cluster_res <- kmeans(cluster_mat, centers = K)

  cluster_data <- data_frame(
    "interp_time" = cluster_data$interp_time,
    "rsv" = cluster_data$rsv,
    "cluster" = cluster_res$cluster
  ) %>%
  mutate_at(vars(-interp_time), as.character) %>%
  mutate_if(is.character, as_factor)

  list("partition" = cluster_data, "cluster" = cluster_res)
}

## ---- data ----
abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 50}, prune = TRUE) %>%
  transform_sample_counts(asinh)

time_mapping <- seq(1, n, length.out = 256)
dup_sample_data <- sample_data(abt)[round(time_mapping), ]
dup_sample_data$interp_time <- 1:256

## ---- dwt ----
x <- t(get_taxa(abt))
n <- nsamples(abt)

detail_coefs <- detail_array(x, scale_lev = TRUE)
cluster_res <- cluster_detail_array(detail_coefs$D, K = 4)

aligned_partition <- cluster_res$partition %>%
  left_join(dup_sample_data) %>%
  arrange(ind, rsv, time)

p <- list()
for (cur_ind in c("F", "D", "E")) {
  p[[cur_ind]] <- ggplot(aligned_partition %>% filter(ind == cur_ind)) +
    geom_tile(aes(x = rsv, y = -interp_time, fill = cluster)) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(condition ~ ., scales = "free") +
    theme(
      panel.spacing = unit(0, "lines"),
      axis.text.x = element_blank(),
      panel.border = element_rect(fill = "transparent", size = 0.5)
    ) +
    ggtitle(sprintf("Subject %s", cur_ind))
}

cluster_res$cluster$centers

## ---- denoising ----
thresh_xwd <- threshold(detail_coefs$transforms[[11]], policy = "cv")
plot(interpolate_dyadic(x[, 11]))
points(wr(thresh_xwd), col = "red")
