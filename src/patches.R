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

ggplot(x[1:100, ]) +
  geom_point(aes(x = time, y = value, group = rsv)) +
  facet_grid(ind ~ .)

times_concat <- x %>%
  group_by(ind) %>%
  select(ind, time) %>%
  unique() %>%
  arrange(ind, time)

times_concat$concat_time <- times_concat$time
for (i in seq_len(nrow(times_concat) - 1)) {
  if (times_concat$concat_time[i + 1] < times_concat$concat_time[i]) {
    cur_diff <- times_concat$concat_time[i] - times_concat$concat_time[i + 1] + 1
    times_concat$concat_time[(i + 1):nrow(times_concat)] <- cur_diff +
      times_concat$concat_time[(i + 1):nrow(times_concat)]
  }
}

x <- x %>%
  left_join(times_concat) %>%
  rename(
    original_time = time,
    time = concat_time
  )

split_rsvs <- dlply(x, .(rsv))
eval_times <- seq(min(x$time), max(x$time), by = 1)
x_interp <- matrix(0, nrow = length(split_rsvs), ncol = length(eval_times))
rownames(x_interp) <- names(split_rsvs)
for (i in seq_along(split_rsvs)) {
  x_interp[i, ] <- approxfun(split_rsvs[[i]]$time, split_rsvs[[i]]$value)(eval_times)
}

x_interp <- x_interp %>%
  matrix_to_df("rsv") %>%
  gather(time, value, -rsv) %>%
  mutate(time = as.integer(gsub("X", "", time))) %>%
  x %>%
  select(sample, ind, time, rsv) %>%
  right_join(x_interp) %>%
  left_join(taxa)

x_mat <- x_interp %>%
  select(time, rsv, value) %>%
  spread(rsv, value) %>%
  select(-time) %>%
  as.matrix()

## ---- define-patches ----
patch_length <- 10
patch_times <- list()
times <- unique(x$time)

for (i in seq_len(length(times) - patch_length)) {
  patch_times[[i]] <- times[i:(i + patch_length - 1)]
}

patch_x <- list()
for (i in seq_along(patch_times)) {
  if (i %% 10 == 0) {
    cat(sprintf("Generating patch %s / %s\n", i, length(patch_times)))
  }
  cur_times <- patch_times[[i]]
  patch_x[[i]] <- x %>%
    filter(time %in% cur_times) %>%
    group_by(ind, rsv) %>%
    summarise(
      patch_start = cur_times[1],
      patch_end = cur_times[length(cur_times)],
      start_condition = condition[1],
      end_condition = condition[length(condition)],
      values = list(setNames(value, seq_along(value)))
    )
}

patch_x <- x %>%
  filter(time %in% cur_times) %>%
  group_by(ind, rsv) %>%
  nest(value)

patch_x <- do.call(rbind, patch_x) %>%
  rownames_to_column("patch_id") %>%
  mutate(
    patch_id = as.integer(patch_id)
  )

patch_df <- unnest(patch_x, .sep = "_")
patch_df$patch_ix <- rep(seq_len(patch_length), length.out = nrow(patch_df))

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
