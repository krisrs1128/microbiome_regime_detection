#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some experiments using wavelet denoising on the antibiotics data, in an
## attempt to cluster the timepoints across different bacteria.

## ---- libraries ----
library("phyloseq")
library("wavethresh")
library("reshape2")
library("tidyverse")
library("forcats")
library("dendextend")
source("data.R")
theme_set(ggscaffold::min_theme())

## ---- utils ----
save_fig <- function(fname, p, output_dir = "../../doc/figure/") {
  ggsave(file.path(output_dir, fname), p)
}

soft_thresh <- function(x, lambda) {
  x[abs(x) < lambda] <- 0
  x[x > lambda] <- x - lambda
  x[x < lambda] <- x + lambda
  x
}

interpolate_dyadic <- function(y) {
  n <- length(y)
  n_interp <- 2 ^ ceiling(log(n, 2))
  f_y <- approxfun(seq_along(y), y, yleft = 0, yright = 0)
  f_y(seq(1, n, length.out = n_interp))
}

stacked_detail <- function(ywd, min_lev = 0, max_lev = NULL, scale_lev = FALSE) {
  if (is.null(max_lev)) {
    max_lev <- ywd$nlevel 
  }

  n_interp <- length(ywd$D) + 1
  levs <- seq(min_lev, max_lev)
  stacked <- list()
  for (l in seq_along(levs)) {
    stacked[[l]] <- rep(accessC(ywd, levs[l]), each = n_interp / 2 ^ l)
    if (levs[l] != max_lev) {
      stacked[[l]] <- cbind(
        stacked[[l]],
        rep(accessD(ywd, levs[l]), each = n_interp / 2 ^ l)
      )
    }
  }

  stacked <- do.call(cbind, stacked)
  if (scale_lev) {
    stacked <- scale(stacked)
    stacked[is.na(stacked)] <- 0
  }

  stacked
}

concat_coef <- function(ywd, min_lev = 0, max_lev = NULL, scale_lev = FALSE) {
  if (is.null(max_lev)) {
    max_lev <- ywd$nlevel
  }
  concat <- list()
  levs <- seq(min_lev, max_lev)
  for (l in seq_along(levs)) {
    if (levs[l] == max_lev) {
      concat[[l]] <- unlist(list("C" = accessC(ywd, levs[l])))
    } else {
      concat[[l]] <- unlist(list("D" = accessD(ywd, levs[l]), "C" = accessC(ywd, levs[l])))
    }

    if (scale_lev) {
      concat[[l]] <- as.numeric(scale(concat[[l]]))
    }
    names(concat[[l]]) <- paste0("l", levs[l], "_", names(concat[[l]]))
  }

  do.call(c, concat)
}

wavelet_transforms <- function(y, ...) {
  p <- ncol(y)
  transforms <- list()
  for (j in seq_len(p)) {
    transforms[[j]] <- x[, j] %>%
      interpolate_dyadic %>%
      wd(...)
  }

  transforms
}

detail_array <- function(transforms, ...) {
  n <- length(transforms[[1]]$D)
  p <- length(transforms)
  D <- array(0, dim = c(2 ^ ceiling(log(n, 2)), p, 2 * ceiling(log(n, 2))))
  for (j in seq_len(p)) {
    D[, j,] <- stacked_detail(transforms[[j]], ...)
  }
  D
}

concat_matrix <- function(transforms, ...) {
  p <- length(transforms)
  C <- list()
  for (j in seq_len(p)) {
    C[[j]] <- concat_coef(transforms[[j]], ...)
  }
 do.call(rbind, C)
}

cluster_detail_array <- function(D, K = 8, scale_clust = FALSE) {
  mdetail <- D %>%
    melt(varnames = c("interp_time", "rsv", "coef_level"))

  cluster_data <- mdetail %>%
    mutate(coef_level = paste0("L", coef_level)) %>%
    spread(coef_level, value)

  cluster_mat <- cluster_data %>%
    select(starts_with("L")) %>%
    as.matrix

  if (scale_clust) {
    cluster_mat <- scale(cluster_mat)
  }
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

replace_coefs <- function(ywd, w, L = 8) {
  for (l in 0:L) {
    lev_label <- sprintf("l%s_", l)
    D_coefs <- w[grep(paste0(lev_label, "D"), names(w))]
    C_coefs <- w[grep(paste0(lev_label, "C"), names(w))]

    coef_order <- order(as.numeric(str_extract(names(C_coefs), "[0-9]+$")))
    if (l < 8) ywd <- putD(ywd, l, D_coefs[coef_order])
    ywd <- putC(ywd, l, C_coefs[coef_order])
  }
  ywd
}


## ---- data ----
abt <- get(load("../../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 50}, prune = TRUE) %>%
  transform_sample_counts(asinh)
x <- t(get_taxa(abt))
n <- nsamples(abt)

time_mapping <- seq(1, n, length.out = 256)
dup_sample_data <- sample_data(abt)[round(time_mapping), ]
dup_sample_data$interp_time <- 1:256

## ---- dwt ----
xwd <- wavelet_transforms(x, family = "DaubExPhase", filter.number = 1)
D <- detail_array(xwd)
cluster_res <- cluster_detail_array(D, K = 4, scale_clust = TRUE)

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
for (i in seq_along(p)) {
  save_fig(sprintf("wavelet_time_cluster-%d.png", i), p[[i]])
}

## ---- centroids ----
z <- x
dimnames(z) <- NULL
mx <- z[time_mapping, ] %>%
  melt(varnames = c("interp_time", "rsv")) %>%
  mutate_at(vars(c("rsv")), funs(as_factor(as.character(.)))) %>%
  left_join(aligned_partition) %>%
  as_data_frame

centroids <- mx %>%
  group_by(interp_time, cluster) %>%
  summarise(
    condition = condition[1],
    mean = mean(value),
    sd = sd(value),
    n_samp = n()
  )

p <- ggplot(centroids) +
  geom_line(aes(x = interp_time, y = mean, col = cluster)) +
  geom_ribbon(aes(x = interp_time, ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd, fill = cluster), alpha = 0.1) +
  facet_wrap(~cluster)
save_fig("wavelet_centroids-sd.png", p)

p <- ggplot(centroids) +
  geom_line(aes(x = interp_time, y = mean, col = cluster)) +
  geom_ribbon(aes(x = interp_time, ymin = mean - 10 / sqrt(n_samp), ymax = mean + 10 / sqrt(n_samp), fill = cluster), alpha = 0.4) +
  facet_grid(cluster ~ ., scale = "free", space = "free")
save_fig("wavelet_centroids-n.png", p)

## ---- denoising ----
x_thresh <- z[time_mapping, ]
for (j in seq_len(ncol(x_thresh))) {
  x_thresh[, j] <- wr(threshold(xwd[[j]], policy = "cv"))
}

## ---- cluster ----
samples <- sample_data(abt)[time_mapping] %>%
  data.frame() %>%
  rownames_to_column("sample")

taxa <- abt %>%
  tax_table %>%
  taxa_labels

colnames(x_thresh) <- colnames(x)
rownames(x_thresh) <- samples$sample
hclust_x <- hclust(dist(t(x_thresh)))

dendro <- as.dendrogram(hclust_x)
dendro <- reorder(dendro, -colMeans(x_thresh))
joined_data <- join_sources(x_thresh, taxa, samples, dendro)

p <- ggplot(joined_data) +
  geom_tile(aes(x = rsv, y = sample, fill = value)) +
  scale_fill_gradient(low = "white", high = "black") +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(ind ~ ., scales = "free", space = "free")
save_fig("stsacked_wavelet_hclust.png", p)

## ---- ts-features ----
xwd_thresh <- lapply(xwd, threshold, policy = "cv")
C <- concat_matrix(xwd_thresh)
clust_C <- hclust(dist(scale(C)))

K <- 4
M <- data.frame(
  "cluster" = as_factor(as.character(cutree(clust_C, K))),
  C
) %>%
  as_data_frame %>%
  gather(coefficient, value, -cluster) %>%
  group_by(cluster, coefficient) %>%
  summarise(value = mean(value)) %>%
  spread(coefficient, value) %>%
  ungroup

centroids <- matrix(0, nrow = K, ncol = length(xwd[[1]]$D) + 1)
for (k in seq_len(K)) {
  centroids[k, ] <- wr(replace_coefs(xwd[[1]], unlist(M[k, -1])))
}

p <- ggplot(melt(centroids, varnames = c("wv_cluster", "interp_time"))) +
  geom_line(aes(x = interp_time, y = value, col = as.factor(wv_cluster))) +
  facet_wrap(~wv_cluster)
save_fig("concat_wavelet_hclust_reconstruction.png", p)

mx$wv_cluster <- as_factor(as.character(cutree(clust_C, K)))
centroids_wv <- mx %>%
  group_by(interp_time, wv_cluster) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )

p <- ggplot(centroids_wv) +
  geom_line(aes(x = interp_time, y = mean, col = wv_cluster)) +
  geom_ribbon(aes(x = interp_time, ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd, fill = wv_cluster), alpha = 0.1) +
  facet_wrap(~wv_cluster)
save_fig("concat_wavelet_hclust_averages.png", p)
