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
data(abt)

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
## read parameter and state info
###############################################################################
emission <- read_csv(
  file.path(opts$dir, "emission.csv"),
  col_names = c("sample", "iter", "K", "param", "param_ix", "value")
)

dynamics <- read_csv(
  file.path(opts$dir, "dynamics.csv"),
  col_names = c("sample", "iter", "K", "param", "param_ix", "value")
)

stateseq <- read_csv(
  file.path(opts$dir, "stateseq.csv"),
  col_names = c("sample", "iter", seq_len(ncol(x_asinh)))
) %>%
  gather(time, K, -sample, -iter)

###############################################################################
## join parameter and state sequence information
###############################################################################
params <- stateseq %>%
  left_join(emission) %>%
  unite(param_param_ix, param, param_ix) %>%
  spread(param_param_ix, value) %>%
  mutate(
    time = as.integer(time),
  )

pcmat <- params %>%
  select(-sample, -iter, -time, -K) %>%
  as.matrix()

pc_res <- princomp(pcmat)

head(pc_res$scores)

pc_df <- cbind(params, pc_res$scores)
ggplot(pc_df) +
  geom_point(aes(x = C_0, y = sqrt(R_0), col = time))
