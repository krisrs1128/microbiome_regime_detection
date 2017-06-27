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
library("mvtnorm")
library("phyloseq")
source("em_hmm.R")

## ---- simulate ----
# underlying states
z <- cbind(
  c(rep(1, 10), rep(2, 4), rep(1, 6), rep(3, 10)),
  c(rep(2, 6), rep(1, 5), rep(4, 5), rep(1, 4), rep(4, 10)),
  c(rep(4, 3), rep(1, 10), rep(2, 5), rep(4, 5), rep(3, 7))
)[rep(1:30, each = 40), c(1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3)]
image(z)
time_len <- nrow(z)
n <- ncol(z)
K <- length(unique(z))

# parameters per state
theta <- list(
  list("mu" = c(0, 0), "sigma" = diag(2)),
  list("mu" = c(-2, -1), "sigma" = diag(2)),
  list("mu" = c(2, 1), "sigma" = diag(2)),
  list("mu" = c(-1, -1), "sigma" = diag(2))
)

# observed data
y <- array(dim = c(time_len, n, 2))
for (k in seq_along(theta)) {
  y[z == k] <- rmvnorm(
    sum(z == k),
    theta[[k]]$mu,
    theta[[k]]$sigma
  )
}

## ---- visualize-simulation ----
image(z)
image(y[,, 1])
image(y[,, 2])
perm_ix <- sample(n)
image(z[, perm_ix])
image(y[, perm_ix, 1])

## ---- hmm ----
em_res <- hmm_em(y, K = 4, n_iter = 10)
em_res$theta
em_z <- apply(em_res$gamma, c(1, 3), which.max)
table(em_z, z)
image(em_z)
dev.new()
image(z)

plot(y[, 2], col =  z[, 2])
plot(y[, 2], col =  apply(em_res$gamma[,, 2], 1, which.max))
em_res$theta
