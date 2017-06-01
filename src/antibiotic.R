#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some brainstorming with the antibiotics data.

## ---- libraries ----
library("tidyverse")
library("proxy")
library("phyloseq")
library("grid")
library("vegan")
library("dendextend")
theme_set(ggscaffold::min_theme(list(
                        "legend_position" = "right",
                        "border_size" = 0.8
                      ))
          )

## --- utils ----
melted_counts <- function(x) {
  x %>% data.frame() %>%
    rownames_to_column("sample") %>%
    as_data_frame %>%
    gather(key = "rsv", value = "value", -sample) %>%
    mutate(
      scaled = asinh(value) / nrow(x),
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

join_sources <- function(x, taxa, samples, dendro, h = 0.5) {
  leaf_ix <- order.dendrogram(dendro)
  leaf_order <- labels(dendro)[leaf_ix]
  cluster <- data.frame("cluster" = cutree(dendro, h = h)) %>%
    rownames_to_column(var = "rsv")

  cluster <- cluster %>%
    left_join(data_frame(rsv = leaf_order, leaf_ix = leaf_ix)) %>%
    arrange(leaf_ix)
  cluster$cluster <- factor(cluster$cluster, levels = unique(cluster$cluster))

  mx <- melted_counts(x) %>%
    left_join(taxa) %>%
    left_join(samples) %>%
    left_join(cluster)

  mx %>%
    group_by(sample, cluster) %>%
    mutate(
      centroid = median(scaled),
      centroid_prob = mean(value)
    )
}

combined_heatmap <- function(mx, fill_type = "bw") {
  p1 <- ggplot(mx) +
    geom_tile(aes(x = leaf_ix, y = sample, fill = scaled)) +
    facet_grid(ind ~ cluster, scales = "free", space = "free") +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      axis.title = element_blank(),
      axis.text = element_blank()
    )

  if (fill_type == "bw") {
    p1 <- p1 + scale_fill_gradient(low = "white", high = "black")
  } else if (fill_type == "gradient2"){
    p1 <- p1 + scale_fill_gradient2(high = "#48b987", low = "#b94e48")
  }

  unique_mx <- mx %>%
    filter(sample == mx$sample[1])
  rep_ix <- rep(1:10, nrow(unique_mx))
  inv_rep_ix <- rep(seq_len(nrow(unique_mx)), each = 10)
  rsvs <- data_frame(
    "leaf_ix" = unique_mx$leaf_ix[inv_rep_ix],
    "cluster" = unique_mx$cluster[inv_rep_ix],
    "y" = rep_ix,
    "label" = unique_mx$label[inv_rep_ix],
    "dummy" = 1
  )

  p2 <- ggplot(rsvs) +
    geom_tile(aes(x = leaf_ix, y = y, fill = label)) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(nrow = 8)) +
    scale_x_continuous(expand = c(0, 0)) +
    facet_grid(dummy ~ cluster, scales = "free", space = "free") +
    theme(
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      strip.text = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      plot.margin = unit(c(0, 0, 0, 0), "null")
    )

  grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
}

centroid_plot <- function(mx) {
  ggplot(mx) +
    geom_point(
      aes(x = time, y = scaled, col = ind), size = 0.2, alpha = 0.2
    ) +
    geom_point(
      aes(x = time, y = centroid), size = 0.7
    ) +
    geom_point(
      aes(x = time, y = centroid, group = rsv, col = ind), size = 0.4
    ) +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~cluster)
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
combined_heatmap(mx)
centroid_plot(mx)

## ---- heatmap-extremes ----
alpha <- 0
D_mix <- alpha * D_jaccard + (1 - alpha) * D_euclidean
tree <- hclust(D_mix, method = "complete")
dendro <- reorder(as.dendrogram(tree), -colMeans(x))
mx <- join_sources(x, taxa, samples, dendro, h = 0.2)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
combined_heatmap(mx)
centroid_plot(mx)

alpha <- 1
D_mix <- alpha * D_jaccard + (1 - alpha) * D_euclidean
tree <- hclust(D_mix)
dendro <- reorder(as.dendrogram(tree), -colMeans(x))
mx <- join_sources(x, taxa, samples, dendro)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
combined_heatmap(mx)
centroid_plot(mx)

## ---- innovations ----
diff_x <- apply(x_scaled, 2, diff)
tree <- hclust(dist(t(diff_x)))
dendro <- reorder(as.dendrogram(tree), -var(diff_x))
mx <- join_sources(diff_x, taxa, samples, dendro, h = 0.15)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
combined_heatmap(mx, "gradient2") 

## ---- innovations-bin ----
diff_x <- apply(x_bin, 2, diff)
tree <- hclust(dist(t(diff_x), method = "Jaccard"))
dendro <- reorder(as.dendrogram(tree), -var(diff_x))
mx <- join_sources(diff_x, taxa, samples, dendro, h = 0.985)
sort(table(mx$cluster), decreasing = TRUE) / nrow(mx)
combined_heatmap(mx, "gradient2") 

## ---- pacf ----
head(x_scaled)
pacfs <- apply(x_scaled, 1, pacf, plot = FALSE)
