#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Some examples using the different HMMs on toy data
##
## author: sankaran.kris@gmail.com
## date: 07/05/2017

###############################################################################
## setup libraries and themes
###############################################################################
library("mvnfast")
library("jsonlite")
library("MCMCpack")
library("Rcpp")
sourceCpp("messages.cpp")
source("utils.R")
library("tidyverse")
library("viridis")
library("reshape2")
theme_set(ggscaffold::min_theme(list(
                        "legend_position" = "right",
                        "border_size" = 0.2
                      ))
          )

sim <- simulate_data()
my <- melt(
  sim$y,
  varnames = c("time", "sample", "p")
)
mz <- melt(
  sim$z,
  varnames = c("time", "sample")
)

p <- ggplot(mz) +
  geom_tile(
    aes(x = time, y = sample, fill = value)
  ) +
  scale_fill_viridis(option = "magma") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none")
ggsave("~/Desktop/lab_meetings/20170705/figure/sticky_hmm_truth.pdf", p, height = 2, width = 4)

p <- ggplot(my %>% filter(p == 1)) +
  geom_tile(
    aes(x = time, y = sample, fill = value)
  ) +
  scale_fill_viridis(option = "magma") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none")
ggsave("~/Desktop/lab_meetings/20170705/figure/sticky_hmm_observed.pdf", p, height = 2, width = 4)

source("bayes_hmm2.R")
res <- block_sampler(sim$y, list(K = 4, n_iter = 100, kappa = 0))
mu1 <- sapply(res$theta, function(x) {x$mu[1]})
mz_est <- melt(
  res$z,
  varnames = c("time", "sample")
)
mz_est$value <- factor(mz_est$value, order(mu1))
mz_est$value <- as.numeric(mz_est$value)

p <- ggplot(mz_est) +
  geom_tile(
    aes(x = time, y = sample, fill = value)
  ) +
  scale_fill_viridis(option = "magma") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none")
ggsave("~/Desktop/lab_meetings/20170705/figure/sticky_hmm_no_kappa.pdf", p, height = 2, width = 4)

## larger kappa now
res <- block_sampler(sim$y, list(K = 4, n_iter = 100, kappa = 3))
mu1 <- sapply(res$theta, function(x) {x$mu[1]})
mz_est <- melt(
  res$z,
  varnames = c("time", "sample")
)
mz_est$value <- factor(mz_est$value, order(mu1))
mz_est$value <- as.numeric(mz_est$value)

p <- ggplot(mz_est) +
  geom_tile(
    aes(x = time, y = sample, fill = value)
  ) +
  scale_fill_viridis(option = "magma") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none")
ggsave("~/Desktop/lab_meetings/20170705/figure/sticky_hmm_large_kappa.pdf", p, height = 2, width = 4)

source("hdp_hmm2.R")
hyper <- list(
  n_iter = 200,
  L = 10,
  gamma = 1e-4, ## controls total number of states
  alpha = 1e-2, ## controls number of states any one state can transition to
  kappa = 1e-2
)
res <- block_sampler(sim$y, hyper)

mu1 <- sapply(res$theta, function(x) {x$mu[1]})
mz_est <- melt(
  res$z,
  varnames = c("time", "sample")
)
mz_est$value <- factor(mz_est$value, order(mu1))
mz_est$value <- as.numeric(mz_est$value)

p <- ggplot(mz_est) +
  geom_tile(
    aes(x = time, y = sample, fill = value)
  ) +
  scale_fill_viridis(option = "magma")  +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = "none")

ggsave("~/Desktop/lab_meetings/20170705/figure/sticky_hdp_hmm.pdf", p, height = 2, width = 4)

