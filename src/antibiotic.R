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
library("depmixS4")
library("abind")
theme_set(ggscaffold::min_theme(list(
                        "legend_position" = "right",
                        "border_size" = 0.2
                      ))
          )
figure_dir <- file.path("..", "doc", "figure")
dir.create(figure_dir)

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
    mutate(centroid_prob = mean(present)) %>%
    group_by(sample, cluster, present) %>%
    mutate(centroid = mean(scaled))
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

  arrangeGrob(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
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
    facet_wrap(~cluster)

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
    facet_wrap(~cluster)

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

## ---- hmm ----
i <- 100
df <- data.frame(y = x_scaled[, i] + runif(nrow(x_scaled), 0, 0.001))
msp <- depmix(y ~ 1, nstates = 4, data = df, emcontrol = em.control(classification = "soft"))
fm <- fit(msp)

plot_hmm <- cbind(df, fm@posterior) %>%
  rownames_to_column("sample") %>%
  left_join(samples) %>%
  gather(state_prob, prob, starts_with("S", ignore.case = FALSE))

all_plots[["hmm-example"]] <- ggplot(plot_hmm) +
  geom_line(aes(x = time, y = y, col = ind), size = 0.4, alpha = 0.3) +
  geom_point(aes(x = time, y = y, col = ind, size = prob)) +
  facet_grid(state_prob ~ .) +
  scale_size_continuous(range = c(0.05, 2)) +
  scale_color_brewer(palette = "Set1")

## ---- parallel-hmm ----
plot_hmm <- list()
K <- 4

## Loop over every RSV and fit an HMM
for (i in seq_len(ncol(x_scaled))) {
  if (i %% 10 == 0) {
    cat(sprintf("HMM on rsv %s\n", i))
  }

  df <- data_frame(
    "sample" = rownames(x_scaled),
    "y" = x_scaled[, i] + runif(nrow(x_scaled), 0, 0.005),
    "rsv" = colnames(x_scaled)[i]
  )

  msp <- depmix(
    y ~ 1,
   nstates = K,
    data = df,
    emcontrol = em.control(classification = "soft")
  )
  fm <- try(fit(msp, verbose = FALSE))

  if (class(fm) != "try-error") {
    ## try to align states across rsvs, naive approach just sorts by emission mean
    state_order <- c(1, 1 + order(summary(fm)[, 1], decreasing = TRUE))
    ordered_posterior <- fm@posterior[, state_order]
    colnames(ordered_posterior) <- c("state", paste0("S", seq_len(K)))

    plot_hmm[[i]] <- cbind(df, ordered_posterior) %>%
      gather(state_prob, prob, starts_with("S", ignore.case = FALSE)) %>%
      as_data_frame
  }
}

plot_hmm <- do.call(rbind, plot_hmm) %>%
  left_join(samples)
plot_hmm$rsv <- factor(plot_hmm$rsv, levels = mix_tree$label)

for (state in paste0("S", 1:K)) {
  all_plots[[sprintf("hmm-%s", state)]] <- ggplot(plot_hmm %>% filter_(sprintf("state_prob == '%s'", state))) +
    geom_tile(aes(x = rsv, y = sample, fill = prob)) +
    scale_fill_gradient(low = "white", high = "black") +
    facet_grid(ind ~ ., scale = "free_y") +
    theme(axis.text = element_blank()) +
    ggtitle(sprintf("State %s", state))
}

## write all the plots to file
for (i in seq_along(all_plots)) {
  cur_name <- paste0(names(all_plots)[i], ".png")
  height <- 2.5
  if (grepl("heatmap", cur_name)) {
    height <- 1.5
  } else if (grepl("centroid", cur_name)) {
    height <- 4
  }

  ggsave(
    file.path(figure_dir, cur_name),
    all_plots[[i]],
    dpi = 450,
    height = height, width = 6.5
  )
}

## ---- write-js ---
js_data <- mx %>%
  ungroup() %>%
  arrange(leaf_ix) %>%
  dplyr::select(sample, rsv, label, scaled) %>%
  rename(
    column = rsv,
    row = sample,
    value = scaled
  )

cat(sprintf("var data = %s;", jsonlite::toJSON(js_data)), file = "~/Desktop/100_days/june2/data.js")

phy <- as.phylo(mix_dendro)
phy_df <- data_frame(
  parent = phy$edge[, 1],
  child = phy$edge[, 2],
  edge_length = phy$edge.length
)
phy_df <- rbind(
  phy_df,
  data_frame(parent = "", child = "621", edge_length = 0.1)
)
cat(sprintf("var tree = %s;", jsonlite::toJSON(phy_df)), file = "~/Desktop/100_days/june2/tree.js")
