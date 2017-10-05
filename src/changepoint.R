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
data(abt)

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

###############################################################################
## Filter data to use for modeling and write to file
###############################################################################
abt <- abt %>%
  filter_taxa(function(x) mean(x > 0) > 0.8, prune = TRUE) %>%
  subset_samples(ind == "F")

x_bern <- 1 * t(get_taxa(abt) > 0)
x_asinh <- t(asinh(get_taxa(abt)))

write_csv(
  data.frame(x_bern),
  file.path("data", "abt_bern.csv"),
  col_names = FALSE
)
write_csv(
  data.frame(x_asinh),
  file.path("data", "abt_asinh.csv"),
  col_names = FALSE
)

###############################################################################
## read samples from different methods
###############################################################################
samples <- read_csv(
  file.path("data", "changepoint", "samples.csv"),
  col_names = c("iter", "changepoint", "row")
)

samples_stats <- samples %>%
  group_by(changepoint, row) %>%
  summarise(n_draws = n())

ggplot(samples_stats) +
  geom_point(
    aes(x = changepoint, y = row, alpha = n_draws)
  ) +
  scale_alpha(range = c(0.05, 1))
