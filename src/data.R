
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
  taxa$label[is.na(taxa$label)] <- "other"
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
    mutate(centroid_prob = mean(present)) %>%
    group_by(sample, cluster, present) %>%
    mutate(centroid = mean(scaled))
}
