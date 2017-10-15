#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Example using SLDS for a single series.
##
## author: sankaran.kris@gmail.com
## date: 10/13/2017

library("tidyverse")
library("phyloseq")
library("treelapse")
library("forcats")
data(abt)

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2", na.value = "black")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2", na.value = "black")

theme_set(theme_bw())
theme_update(
  panel.border = element_rect(size = 0.5),
  panel.background = element_rect(fill = "#F7F7F7"),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

opts <- list(
  "dir" = file.path("..", "..", "data", "slds"),
  "k_filter" = 0.2
)
dir.create(opts$dir)

###############################################################################
## Filter data to use for modeling and write to file
###############################################################################
data(abt)
abt <- abt %>%
  subset_samples(ind == "F") %>%
  filter_taxa(function(x) mean(x > 0) > opts$k_filter, prune = TRUE)

x_asinh <- asinh(get_taxa(abt))
write_csv(
  data.frame(x_asinh),
  file.path(opts$dir, "abt.csv"),
  col_names = FALSE
)

###############################################################################
## taxonomic and sample information
###############################################################################
taxa <- tax_table(abt) %>%
  data.frame() %>%
  rownames_to_column("seq")
#taxa$family <- substr(taxa$Taxon_5, 0, 8)
taxa$family <- taxa$Taxon_5
taxa$family[taxa$family == ""] <- NA
taxa$family <- fct_lump(taxa$family, 7)
taxa$family <- taxa$family %>%
  recode(
    Alcaligenaceae_Sutterella = "Sutterella",
    Peptostreptococcaceae_1 = "Peptostreptococcaceae"
  )
taxa$family[is.na(taxa$family)] <- "Other"

sample_df <- sample_data(abt) %>%
  data.frame() %>%
  rownames_to_column("sample")

###############################################################################
## read parameter and state info
###############################################################################
emission <- read_csv(
  file.path(opts$dir, "emission.csv"),
  col_names = c("seq_ix", "iter", "K", "param", "param_ix", "value")
)

dynamics <- read_csv(
  file.path(opts$dir, "dynamics.csv"),
  col_names = c("seq_ix", "iter", "K", "param", "param_ix", "value")
)

stateseq <- read_csv(
  file.path(opts$dir, "stateseq.csv"),
  col_names = c("seq_ix", "iter", seq_len(ncol(x_asinh)))
) %>%
  gather(sample_ix, K, -seq_ix, -iter)

emit_dyn <- emission %>%
  bind_rows(dynamics)

###############################################################################
## join parameter and state sequence information
###############################################################################
params <- stateseq %>%
  left_join(emit_dyn) %>%
  unite(param_param_ix, param, param_ix) %>%
  spread(param_param_ix, value) %>%
  mutate(
    sample_ix = as.integer(sample_ix),
    R_0 = sqrt(R_0),
    Q_0 = sqrt(Q_0)
  ) %>%
  mutate(
    sample = colnames(x_asinh)[sample_ix],
    seq = rownames(x_asinh)[1 + seq_ix]
  ) %>%
  left_join(taxa) %>%
  left_join(sample_df) %>%
  mutate(
    time_bin = cut(time, seq(0, 60, by = 10), include.lowest = TRUE)
  )

###############################################################################
## Plot scores and loadings from PC on this matrix
###############################################################################
pcmat <- params %>%
  select(
    starts_with("A_"),
    starts_with("C_"),
    starts_with("Q_"),
    starts_with("R_")
  ) %>%
  as.matrix()

pcmat[is.na(pcmat)] <- 0

pc_res <- princomp(scale(pcmat))
pc_df <- cbind(params, pc_res$scores)

ggplot(pc_df %>%
       filter(iter > 150)
       ) +
  geom_vline(xintercept = 0, alpha = 0.4, size = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5) +
  geom_point(
    aes(x = Comp.1, y = Comp.2, col = family),
    alpha = 0.1, size = 0.1
  ) +
  facet_grid(time_bin ~ family) +
  scale_color_brewer(
    palette = "Set3",
    guide = guide_legend(override.aes = list("alpha" = 1, "size" = 1))
  ) +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_text(angle = 0),
    panel.spacing.x = unit(0, "cm"),
    legend.position = "bottom"
  ) +
  coord_fixed(sqrt(pc_res$sdev[2] / pc_res$sdev[1])) +
  xlim(-3.5, 4.2) +
  ylim(-3, 4)

ggsave("../../doc/figure/slds_pca_scores.png", width = 6.5, height = 4.9)

## plot the loadings
pc_load <- loadings(pc_res)
class(pc_load) <- "matrix"
pc_load <- melt(pc_load) %>%
  spread(Var2, value)

ggplot(pc_load) +
  geom_text(
    aes(x = Comp.1, y = Comp.2, size = Comp.3, label = Var1)
  ) +
  geom_vline(xintercept = 0, alpha = 0.4, size = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.5) +
  geom_segment(
    aes(x = 0, y = 0, xend = 0.9 * Comp.1, yend = 0.9 * Comp.2),
    size = 0.5, alpha = 0.9
  ) +
  scale_size(range = c(2, 5)) +
  coord_fixed(sqrt(pc_res$sdev[2] / pc_res$sdev[1])) +
  theme(
    legend.position = "bottom"
  )
ggsave("../../doc/figure/slds_pca_loadings.png")

###############################################################################
## heatmap of posterior means
###############################################################################
params$NA_NA <- NULL
hm_df <- params %>%
  gather(param, value, A_0, C_0, Q_0, R_0) %>%
  group_by(param) %>%
  mutate(value = value / max(value, na.rm = TRUE)) %>%
  group_by(seq, time, param) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    family = family[1]
  )

seq_order <- pc_df %>%
  group_by(seq) %>%
  summarise(c1 = mean(Comp.1, na.rm = TRUE)) %>%
  arrange(c1) %>%
  .[["seq"]]

hm_df$seq <- factor(
  hm_df$seq,
  levels = seq_order
)

## threshold some outliers
hm_df$value[hm_df$value < -0.08] <- -0.08
hm_df$value[hm_df$value > 0.08] <- 0.08

hm_df <- hm_df %>%
  ungroup() %>%
  mutate(
    family = recode(
      family,
      Alistipes = "Al.",
      Ruminococcaceae = "Rumino.",
      Lachnospiraceae = "Lachno.",
      Streptococcaceae = "Strep.",
      Parabacteroides = "Parab.",
      Peptostreptococcaceae = "PStrep.",
      Suterella = "Sutter.",
      Veillonellaceae = "Veill."
    )
  )

ggplot(hm_df) +
  geom_tile(
    aes(x = seq, y = time, fill = value)
  ) +
  facet_grid(param ~ family, scale = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 20, 40)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#78749a", mid = "#f7f7f7", high = "#59ac53",
    guide = guide_colorbar(barheight = 0.3, ticks = FALSE)
  ) +
  theme(
    panel.spacing.x = unit(0, "cm"),
    axis.text.x = element_blank(),
    strip.text.x = element_text(angle = 90, hjust = 0, size = 5),
    legend.position = "bottom"
  )

ggsave(
  "../../doc/figure/slds_parameter_heatmap.png",
  width = 8.42,
  height = 3.08
)
