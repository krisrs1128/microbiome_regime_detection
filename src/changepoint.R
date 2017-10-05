#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Changepoint detection using Zhou and Lester's BASIC model.
##
## author: sankaran.kris@gmail.com
## date: 10/05/2017

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("phyloseq")
library("treelapse")
library("forcats")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_update(
  panel.border = element_rect(size = 0.5),
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
  "dir" = file.path("..", "data", "changepoint"),
  "k_filter" = 0.07
)

###############################################################################
## Filter data to use for modeling and write to file
###############################################################################
data(abt)
abt <- abt %>%
  subset_samples(ind == "F") %>%
  filter_taxa(function(x) mean(x > 0) > opts$k_filter, prune = TRUE) 

x_bern <- 1 * t(get_taxa(abt) > 0)
x_asinh <- asinh(get_taxa(abt))

write_csv(
  data.frame(x_bern),
  file.path(opts$dir, "abt_bern.csv"),
  col_names = FALSE
)
write_csv(
  data.frame(x_asinh),
  file.path(opts$dir, "abt_asinh.csv"),
  col_names = FALSE
)

system("python changepoint.py")

###############################################################################
## Supplementary taxa data
###############################################################################
taxa <- tax_table(abt) %>%
  data.frame() %>%
  rownames_to_column("seq") %>%
  mutate_all(.funs = list(as.character)) %>%
  mutate(family = Taxon_5)

taxa$family[taxa$family == ""] <- taxa$Taxon_4[taxa$family == ""]
taxa$family <- as.factor(taxa$family)
taxa$family <- fct_lump(taxa$family, 7)

taxa$family <- factor(
  taxa$family,
  levels = names(sort(table(taxa$family), decreasing = TRUE))
)


###############################################################################
## read samples from different methods
###############################################################################
samples <- read_csv(
  file.path(opts$dir, "samples.csv"),
  col_names = c("iter", "changepoint", "row")
)

samples$seq <- taxa_names(abt)[1 + samples$row]
samples$time <- sample_data(abt)$time[1 + samples$changepoint]

## sort sequences according to a clustering on the changepoint locations
change_indic <- samples %>%
  dcast(seq ~ time)
D <- dist(log(1 + change_indic[, -1]))
samples$seq <- factor(
  samples$seq,
  levels = change_indic$seq[hclust(D)$order]
)
taxa$seq <- factor(taxa$seq, levels = levels(samples$seq))

sample_stats <- samples %>%
  group_by(seq, time) %>%
  summarise(n_draws = n()) %>%
  left_join(taxa)

ggplot(sample_stats) +
  geom_tile(
    aes(x = time, y = seq, alpha = n_draws, fill = family)
  )

q_values <- read_csv(
  file.path("data", "")
)
