#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some brainstorming with the antibiotics data.

## ---- libraries ----
library("tidyverse")
library("proxy")
library("phyloseq")
library("grid")
library("gridExtra")
library("vegan")
library("ape")
library("dendextend")
library("viridis")
library("ggscaffold")
source("data.R")
theme_set(min_theme(
  list(
    "legend_position" = "right",
    "border_size" = 0.2
  )
))
figure_dir <- file.path("..", "doc", "figure")
dir.create(figure_dir)

## --- utils ----
combined_heatmap <- function(mx, fill_type = "bw") {
  p1 <- ggplot(mx) +
    geom_tile(aes(y = sample, x = leaf_ix, fill = scaled)) +
    facet_grid(ind ~ cluster, scales = "free", space = "free") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      strip.text.x = element_text(size = 4),
    )

  if (fill_type == "bw") {
    p1 <- p1 + scale_fill_viridis(option = "magma", direction = -1)
  } else if (fill_type == "gradient2"){
    p1 <- p1 + scale_fill_gradient2(high = "#32835f", low = "#833256")
  }

  unique_mx <- mx %>%
    filter(sample == mx$sample[1])
  rep_ix <- rep(1:10, nrow(unique_mx))
  inv_rep_ix <- rep(seq_len(nrow(unique_mx)), each = 10)
  rsvs <- data_frame(
    "leaf_ix" = unique_mx$leaf_ix[inv_rep_ix],
    "cluster" = unique_mx$cluster[inv_rep_ix],
    "y" = rep_ix,
    "label" = unique_mx$label[inv_rep_ix]
  )

  p2 <- ggplot(rsvs) +
    geom_tile(aes(x = leaf_ix, y = y, fill = label)) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(ncol = 8)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(. ~ cluster, scales = "free", space = "free") +
    theme(
      panel.border = element_blank(),
      legend.position = "bottom",
      axis.title = element_blank(),
      axis.text = element_blank(),
      strip.text = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      plot.margin = unit(c(0, 0, 0, 0), "null")
    )

  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  rm_ix <- g1$layout$r[grep("strip-r", g1$layout$name)]
  arrangeGrob(rbind(g1[, -rm_ix], g2))
}

centroid_plot <- function(mx) {
  centroid_data <- mx %>%
    group_by(sample, cluster) %>%
    summarise(
      time = time[1],
      ind = ind[1],
      centroid = centroid[1],
      centroid_prob = centroid_prob[1]
    )

  p1 <- ggplot(mx) +
    geom_point(
      aes(x = time, y = scaled, col = ind), size = 0.2, alpha = 0.2
    ) +
    geom_point(
      data = centroid_data,
      aes(x = time, y = centroid), size = 0.7
    ) +
    geom_point(
      data = centroid_data,
      aes(x = time, y = centroid, col = ind), size = 0.4
    ) +
    scale_color_brewer(palette = "Set1") +
    theme(axis.text.x = element_blank()) +
    facet_wrap(~cluster, ncol = 10)

  p2 <- ggplot(mx) +
    geom_point(
      aes(x = time, y = present, col = ind),
      size = 0.2, alpha = 0.2,
      position = position_jitter(height = 0.2)
    ) +
    geom_point(
      data = centroid_data,
      aes(x = time, y = centroid_prob), size = 0.7
    ) +
    geom_point(
      data = centroid_data,
      aes(x = time, y = centroid_prob, col = ind), size = 0.4
    ) +
    scale_color_brewer(palette = "Set1") +
    theme(axis.text.x = element_blank()) +
    facet_wrap(~cluster, ncol = 10)

  list(
    "conditional" = p1,
    "presence" = p2
  )
}

## ---- data ----
download.file("https://github.com/krisrs1128/treelapse/raw/master/data/abt.rda", "../data/abt.rda")
abt <- get(load("../data/abt.rda")) %>%
  filter_taxa(function(x) { var(x) > 5 }, TRUE)

x <- t(get_taxa(abt))

## ---- pos-bin ----
## separate into conditional positive and present absence data
x_bin <- x > 0
class(x_bin) <- "numeric"

x_pos <- x
x_pos[x_pos == 0] <- NA

## ---- distances ----
D_jaccard <- dist(t(x_bin), method = "jaccard")
jaccard_tree <- hclust(D_jaccard)
plot(jaccard_tree)

x_scaled <- asinh(x) / nrow(x)

D_euclidean <- dist(t(x_scaled), method = "euclidean")
euclidean_tree <- hclust(D_euclidean)
plot(euclidean_tree)

alpha <- 0.5
D_mix <- alpha * D_jaccard + (1 - alpha) * D_euclidean
mix_tree <- hclust(D_mix)
plot(mix_tree)

## ---- heatmap-mix ----
samples <- sample_data(abt) %>%
  data.frame() %>%
  rownames_to_column("sample")
taxa <- abt %>%
  tax_table %>%
  taxa_labels

mix_dendro <- reorder(as.dendrogram(mix_tree), -colMeans(x))
mx <- join_sources(x, taxa, samples, mix_dendro, h = 0.5)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)

## save figures into list
all_plots <- list()
all_plots[["heatmap-mix"]] <- combined_heatmap(mx)
p <- centroid_plot(mx)
all_plots[["centroid-mix-conditional"]] <- p$conditional
all_plots[["centroid-mix-presence"]] <- p$presence

## ---- heatmap-extremes ----
alpha <- 0
D_mix <- alpha * D_jaccard + (1 - alpha) * D_euclidean
tree <- hclust(D_mix, method = "complete")
dendro <- reorder(as.dendrogram(tree), -colMeans(x))
mx <- join_sources(x, taxa, samples, dendro, h = 0.2)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
all_plots[["heatmap-euclidean"]] <- combined_heatmap(mx)
p <- centroid_plot(mx)
all_plots[["centroid-euclidean-conditional"]] <- p$conditional
all_plots[["centroid-euclidean-presence"]] <- p$presence

alpha <- 1
D_mix <- alpha * D_jaccard + (1 - alpha) * D_euclidean
tree <- hclust(D_mix)
dendro <- reorder(as.dendrogram(tree), -colMeans(x))
mx <- join_sources(x, taxa, samples, dendro, h = 0.9)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
all_plots[["heatmap-jaccard"]] <- combined_heatmap(mx)
p <- centroid_plot(mx)
all_plots[["centroid-jaccard-conditional"]] <- p$conditional
all_plots[["centroid-jaccard-presence"]] <- p$presence

## ---- innovations ----
diff_x <- apply(x_scaled, 2, diff)
tree <- hclust(dist(t(diff_x)))
dendro <- reorder(as.dendrogram(tree), -var(diff_x))
mx <- join_sources(diff_x, taxa, samples, dendro, h = 0.15)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
all_plots[["heatmap-innovations"]] <- combined_heatmap(mx, "gradient2")
p <- centroid_plot(mx)
all_plots[["centroid-innovations-conditional"]] <- p$conditional
all_plots[["centroid-innovations-presence"]] <- p$presence

## ---- innovations-bin ----
diff_x <- apply(x_bin, 2, diff)
tree <- hclust(dist(t(diff_x), method = "Jaccard"))
dendro <- reorder(as.dendrogram(tree), -var(diff_x))
mx <- join_sources(diff_x, taxa, samples, dendro, h = 0.985)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
all_plots[["heatmap-innovations-bin"]] <- combined_heatmap(mx, "gradient2")
p <- centroid_plot(mx)
all_plots[["centroid-innovations-bin-conditional"]] <- p$conditional
all_plots[["centroid-innovations-bin-presence"]] <- p$presence

## ---- pacf ----
pacfs <- t(apply(x_scaled, 2, function(x) { pacf(x, plot = FALSE)$acf[, 1, 1] }))
pacfs <- pacfs %>%
  data.frame %>%
  rownames_to_column("rsv") %>%
  gather(lag, "value", -rsv) %>%
  as_data_frame() %>%
  mutate(lag = as.integer(gsub("X", "", lag))) %>%
  left_join(taxa)

all_plots[["pacf"]] <- ggplot(pacfs) +
  geom_line(
    aes(x = lag, y = value, group = rsv, col = label),
    alpha = 0.2, size = 0.4
  ) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~label)

## ---- time-pairs ----
x_pairs <- x_scaled
x_offset <- rbind(x_scaled[-1, ], NA)
colnames(x_offset) <- paste0(colnames(x_offset), "_next")
x_pairs <- cbind(x_scaled, x_offset) %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  gather(rsv, "value", -sample) %>%
  mutate(
    next_val = ifelse(grepl("_next", rsv), "post", "pre"),
    rsv = gsub("_next", "", rsv)
  ) %>%
  spread(next_val, value) %>%
  as_data_frame() %>%
  left_join(samples) %>%
  left_join(taxa)

all_plots[["pacf_pairs"]] <- ggplot(x_pairs %>% filter(time >= 10, time <= 20)) +
  geom_abline(slope = 1, size = 1, alpha = 0.3) +
  geom_point(
    aes(x = pre, y = post, col = label),
    alpha = 0.2, size = 0.2
  ) +
  coord_fixed() +
  scale_color_brewer(palette = "Set2") +
  facet_grid(ind~time)

## write all the plots to file
for (i in seq_along(all_plots)) {
  cur_name <- paste0(names(all_plots)[i], ".png")
  width <- 2.5
  if (grepl("heatmap", cur_name)) {
    width <- 3
  } else if (grepl("centroid", cur_name)) {
    width <- 6.5
  }

  ggsave(
    file.path(figure_dir, cur_name),
    all_plots[[i]],
    dpi = 900,
    height = 5.5, width = width
  )
}
