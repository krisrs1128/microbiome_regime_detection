#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some experiments viewing "patches" of time series associated with the
## microbiome.

## ---- libraries ----
library("phyloseq")
library("PMA")
library("tidyverse")
theme_set(ggscaffold::min_theme())

## ---- utils ----
matrix_to_df <- function(x, row_id) {
  x %>%
    data.frame %>%
    rownames_to_column(row_id) %>%
    as_data_frame
}

## ---- get-data ----
download.file("https://github.com/krisrs1128/treelapse/raw/master/data/abt.rda", "../data/abt.rda")
abt <- get(load("../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5 }, TRUE) %>%
  transform_sample_counts(asinh)

x <- t(get_taxa(abt)) %>%
  matrix_to_df("sample") %>%
  gather(rsv, value, -sample)

taxa <- tax_table(abt) %>%
  matrix_to_df("rsv")

samples <- sample_data(abt) %>%
  matrix_to_df("sample")

x <- x %>%
  left_join(taxa) %>%
  left_join(samples)

rsvs <- unique(x$rsv)
inds <- unique(x$ind)

x_interp <- vector(mode = "list", length = length(inds))
names(x_interp) <- inds

for (i in seq_along(inds)) {
  cur_times <- x %>%
    filter(ind == inds[i]) %>%
    .[["time"]] %>%
    unique()
  time_grid <- seq(min(cur_times), max(cur_times), by = 1)
  x_interp[[i]] <- matrix(0, nrow = length(rsvs), ncol = length(time_grid))
  rownames(x_interp[[i]])<- rsvs
  colnames(x_interp[[i]]) <- time_grid

  for (j in seq_along(rsvs)) {
    cat(sprintf("processing rsv %s\n", rsvs[j]))
    cur_data <- x %>%
      filter(ind == inds[i], rsv == rsvs[j]) %>%
      select(time, value)
    x_interp[[i]][j, ] <- approxfun(cur_data$time, cur_data$value)(time_grid)
  }
}

x_interp_df <- x_interp %>%
  matrix_to_df("rsv") %>%
  gather(time, value, -rsv) %>%
  separate(time, c("ind", "time"), sep = "\\.") %>%
  left_join(taxa)

## ---- define-patches ----
patch_length <- 5
patches <- list()

for (i in seq_along(x_interp)) {
  patches[[i]] <- list()
  start_cols <- seq(1, ncol(x_interp[[i]]) - patch_length, by = ceiling(0.5 * patch_length))
  for (j in seq_along(start_cols)) {
    cur <- start_cols[j]
    patches[[i]][[j]] <- x_interp[[i]][, cur:(cur + patch_length - 1)]
    colnames(patches[[i]][[j]]) <- colnames(x_interp[[i]])[cur:(cur + patch_length - 1)]
  }
}

plot(patches[[1]][[1]][4, ])
plot(patches[[1]][[2]][4, ])
dimnames(patches[[1]][[1]])
dimnames(patches[[1]][[2]])

test <- do.call(rbind, patches[[1]])
heatmap(test)

centroids <- kmeans(test, 10)$centers
plot(centroids[1, ])
plot(centroids[2, ])
plot(centroids[3, ])

## ---- pca ----
patch_mat <- do.call(rbind, patch_x$values)
pca_patches <- svd(patch_mat)

W <- pca_patches$u %*% diag(pca_patches$d)
D <- pca_patches$v

D <- D %>%
  matrix_to_df("patch_ix") %>%
  gather(k, value, -patch_ix) %>%
  mutate(
    k = as.numeric(gsub("X", "", k))
  )

ggplot(D) +
  geom_line(aes(x = patch_ix, y = value, group = k)) +
  facet_wrap(~k)

W <- W %>%
  matrix_to_df("patch_id") %>%
  gather(k, weight, -patch_id) %>%
  mutate(
    k = as.integer(gsub("X", "", k)),
    patch_id = as.integer(gsub("X", "", patch_id))
  )

patch_dict <- patch_x %>%
  left_join(W) %>%
  left_join(taxa)

ggplot(patch_dict %>% filter(k != 1)) +
  geom_tile(
    aes(x = rsv, y = -patch_start, fill = weight)
  ) +
  facet_grid(k ~ Taxon_4, scales = "free", space = "free") +
  scale_fill_gradient2(low = "orange", high = "steelblue") +
  theme(
    axis.text.x = element_blank()
  )

## ---- pmd ----
pmd_cv <- PMD.cv(t(scale(t(patch_mat), scale = FALSE)), sumabss = seq(0.5, 1, by = 0.05))
pmd_cv
pmd_res <- PMD(t(scale(t(patch_mat), scale = FALSE)), K = 5, sumabs = 0.8)
pmd_res

V <- pmd_res$v %>%
  matrix_to_df("patch_ix") %>%
  gather(k, value, -patch_ix) %>%
  mutate(
    patch_ix = as.integer(gsub("X", "", patch_ix)),
    k = as.integer(gsub("X", "", k))
  )

ggplot(V) +
  geom_line(
    aes(x = patch_ix, y = value)
  ) +
  facet_wrap(~k)

U <- pmd_res$u %*% diag(pmd_res$d) %>%
  matrix_to_df("patch_id") %>%
  gather(k, weight, -patch_id) %>%
  mutate(
    k = as.integer(gsub("X", "", k)),
    patch_id = as.integer(gsub("X", "", patch_id))
  )

U_dict <- patch_x %>%
  left_join(U) %>%
  left_join(taxa) %>%
  arrange(rsv, patch_start)

ggplot(U_dict) +
  geom_tile(
    aes(x = rsv, y = -patch_start, fill = weight)
  ) +
  facet_grid(k ~ Taxon_4, scales = "free", space = "free") +
  scale_fill_gradient2(low = "orange", high = "steelblue") +
  theme(
    axis.text.x = element_blank()
  )

patch_mat_hat <- rowMeans(patch_mat) %*% t(rep(1, patch_length)) + pmd_res$u %*% diag(pmd_res$d) %*% t(pmd_res$v) %>%
  matrix_to_df("patch_id") %>%
  gather(patch_ix, y_hat, -patch_id) %>%
  mutate(
    patch_ix = as.integer(gsub("X", "", patch_ix)),
    patch_id = as.integer(gsub("X", "", patch_id))
  ) %>%
  left_join(patch_df)

ggplot(patch_mat_hat) +
  geom_point(aes(x = values, y = y_hat))
