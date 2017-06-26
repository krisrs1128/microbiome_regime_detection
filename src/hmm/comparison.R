#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Here, we'll describe a few variants of ordinary HMMs on parallel series. The
## ultimate goal is to partition time points and bacteria in dynamic
## perturbation studies, so that we can see how different bacteria respond
## differently to different perturbations.
##
## author: sankaran.kris@gmail.com
## date: 06/26/2017

## ---- setup ----
library("tidyverse")
library("phyloseq")
source("em_hmm.R")

## ---- simulate ----
# underlying states
z <- cbind(
  c(rep(1, 10), rep(2, 4), rep(1, 6), rep(3, 10)),
  c(rep(2, 6), rep(1, 5), rep(4, 5), rep(1, 4), rep(4, 10)),
  c(rep(4, 3), rep(1, 10), rep(2, 5), rep(4, 5), rep(3, 7))
)[, c(rep(1, 5), rep(2, 4), rep(3, 10))]
image(z)
n <- nrow(z)
p <- ncol(z)
k <- length(unique(z))

# parameters per state
theta <- list(
  list("mu" = 0, "sigma" = 1),
  list("mu" = 2, "sigma" = 1),
  list("mu" = 0.2, "sigma" = 2),
  list("mu" = -0.5, "sigma" = 0.5)
)

# observed data
y <- matrix(0, n, p)
for (k in seq_along(theta)) {
  y[z == k] <- rnorm(
    sum(z == k),
    theta[[k]]$mu,
    theta[[k]]$sigma
  )
}

## ---- visualize-simulation ----
image(z)
image(y)
perm_ix <- sample(p)
image(z[, perm_ix])
image(y[, perm_ix])

## ---- hmm ----
em_res <- hmm_em(y, K = 4)
em_res$theta
em_z <- apply(em_res$gamma, c(1, 3), which.max)
table(em_z, z)
image(em_z)
dev.new()
image(z)
