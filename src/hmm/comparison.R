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
source("sim.R")

## ---- simulate ----
# underlying states
sim <- simulate_data()
y <- sim$y
z <- sim$z

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
