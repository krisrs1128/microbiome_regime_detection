#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some brainstorming with the antibiotics data.

## ---- libraries ----
library("tidyverse")
library("proxy")
library("phyloseq")
library("grid")
library("vegan")
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

join_sources <- function(x, taxa, samples, rsv_order) {
  mx <- melted_counts(x) %>%
    left_join(taxa) %>%
    left_join(samples)
  mx$rsv <- factor(mx$rsv, levels = rsv_order)
  mx %>%
    arrange(rsv)
}

combined_heatmap <- function(mx) {
  p1 <- ggplot(mx) +
    geom_tile(aes(x = rsv, y = sample, fill = scaled)) +
    facet_grid(ind ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "white", high = "#0a0a0a") +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      axis.title = element_blank(),
      axis.text = element_blank()
    )

  unique_mx <- mx %>%
    filter(sample == "D1")
  rep_ix <- rep(1:10, nrow(unique_mx))
  inv_rep_ix <- rep(seq_len(nrow(unique_mx)), each = 10)
  rsvs <- data_frame(
    "rsv" = unique_mx$rsv[inv_rep_ix],
    "y" = rep_ix,
    "label" = unique_mx$label[inv_rep_ix],
    "dummy" = 1
  )

  p2 <- ggplot(rsvs) +
    geom_tile(aes(x = rsv, y = y, fill = label)) +
    scale_fill_brewer(palette = "Set2", guide = guide_legend(nrow = 8)) +
    facet_grid(dummy ~ ., scales = "free", space = "free") +
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

## ---- data ----
download.file("https://github.com/krisrs1128/treelapse/raw/master/data/abt.rda", "../data/abt.rda")
abt <- get(load("../data/abt.rda"))
  ## filter_taxa(function(x) { var(x) > 5 }, TRUE)

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

leaf_ix <- order.dendrogram(reorder(as.dendrogram(mix_tree), -colMeans(x)))
mx <- join_sources(x, taxa, samples, mix_tree$label[leaf_ix])
combined_heatmap(mx)

## ---- heatmap-extremes ----
alpha <- 0
D_mix <- alpha * D_jaccard + (1 - alpha) * D_euclidean
mix_tree <- hclust(D_mix)
leaf_ix <- order.dendrogram(reorder(as.dendrogram(mix_tree), -colMeans(x)))
mx <- join_sources(x, taxa, samples, mix_tree$label[leaf_ix])
combined_heatmap(mx)

alpha <- 1
D_mix <- alpha * D_jaccard + (1 - alpha) * D_euclidean
mix_tree <- hclust(D_mix)
leaf_ix <- order.dendrogram(reorder(as.dendrogram(mix_tree), -colMeans(x_bin)))
mx <- join_sources(x, taxa, samples, mix_tree$label[leaf_ix])
combined_heatmap(mx)
