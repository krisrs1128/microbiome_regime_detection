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
  scale_colour_brewer(..., palette="Set2", na.value = "black")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2", na.value = "black")

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
  "dir" = file.path("..", "..", "data", "changepoint"),
  "k_filter" = 0.2
)

#' Prepare Samples for Plotting
process_samples <- function(samples, abt) {
  samples$seq <- taxa_names(abt)[1 + samples$row]
  samples$time <- sample_data(abt)$time[1 + samples$changepoint]
  samples$seq <- factor(
    samples$seq,
    levels = names(sort(taxa_sums(abt), decreasing = TRUE))
  )
  taxa$seq <- factor(taxa$seq, levels = levels(samples$seq))

  samples %>%
    group_by(seq, time) %>%
    summarise(n_draws = n()) %>%
    left_join(taxa)
}

#' Plot Changepoints
plot_samples <- function(sample_stats) {
  ggplot(sample_stats) +
    geom_tile(
      aes(x = seq, y = time, alpha = n_draws, fill = family)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_alpha(range = c(0, 1)) +
    facet_grid(. ~ family, scale = "free", space = "free") +
    theme(
      strip.text = element_blank(),
      panel.spacing = unit(0, "cm"),
      axis.text.x = element_blank(),
      legend.position = "bottom"
    ) 
}

#' Plot pi and q
plot_pi_q <- function(q_vals, pi_q) {
  q_vals$index <- seq_len(nrow(q_vals)) / nrow(q_vals)
  p <- list()
  p[["q_vals"]] <- ggplot(q_vals) +
    geom_bar(
      aes(x = index, y = q),
      width = 1,
      stat = "identity"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))

  ## plot pi_q
  pi_q$q <- seq_len(nrow(pi_q)) / nrow(q_vals)
  p[["pi_q"]] <- ggplot(pi_q) +
    geom_bar(
      aes(x = q, y = pi_q),
      stat = "identity"
    ) +
    scale_y_continuous(limits = c(0, 0.6), expand = c(0, 0))

  p
}

###############################################################################
## Filter data to use for modeling and write to file
###############################################################################
data(abt)
abt <- abt %>%
  subset_samples(ind == "F") %>%
  filter_taxa(function(x) mean(x > 0) > opts$k_filter, prune = TRUE)

x_bern <- 1 * (get_taxa(abt) > 0)
x_asinh <- asinh(get_taxa(abt))

write_csv(
  data.frame(x_bern),
  file.path(opts$dir, "bern/abt.csv"),
  col_names = FALSE
)
write_csv(
  data.frame(x_asinh),
  file.path(opts$dir, "asinh/abt.csv"),
  col_names = FALSE
)

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
## Gaussian mean / variance changepoint detection
###############################################################################
system("python changepoint.py --dir ../../data/changepoint/asinh/ --iter 2500")
samples <- read_csv(
  file.path(opts$dir, "asinh", "samples.csv"),
  col_names = c("iter", "changepoint", "row")
)
samples %>%
  process_samples(abt) %>%
  plot_samples()
ggsave("../../doc/figure/basic_heatmap.png", width = 6.65, height = 3.15)

## plot q values
q_vals <- read_csv(file.path(opts$dir, "asinh", "q_vals.csv"), col_names = "q")
pi_q <- read_csv(file.path(opts$dir, "asinh", "pi_q.csv"), col_names = "pi_q")
p <- plot_pi_q(q_vals, pi_q)
ggsave("../../doc/figure/changepoint_q.png", p[[1]], width = 6.65, height = 3.15)
ggsave("../../doc/figure/changepoint_eb_prior.png", p[[2]], width = 4, height = 3)

###############################################################################
## Bernoulli changepoint detection
###############################################################################
system("python changepoint.py --dir ../../data/changepoint/bern/ --model bernoulli --iter 2500")
samples <- read_csv(
  file.path(opts$dir, "bern", "samples.csv"),
  col_names = c("iter", "changepoint", "row")
)
samples %>%
  process_samples(abt) %>%
  plot_samples()
ggsave("../../doc/figure/basic_bern_heatmap.png", width = 6.65, height = 3.15)

## plot q values
q_vals <- read_csv(file.path(opts$dir, "bern", "q_vals.csv"), col_names = "q")
pi_q <- read_csv(file.path(opts$dir, "bern", "pi_q.csv"), col_names = "pi_q")
plot_pi_q(q_vals, pi_q)
ggsave("../../doc/figure/changepoint_bern_q.png", p[[1]], width = 4, height = 2.5)
ggsave("../../doc/figure/changepoint_bern_eb_prior.png", p[[2]], width = 4, height = 2.5)
