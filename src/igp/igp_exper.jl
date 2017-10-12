#!/usr/bin/env julia

## File description -------------------------------------------------------------
##
## Experiment using the mixture of GPs sampler on data drawn from the assumed
## model.
##
## author: sankaran.kris@gmail.com
## date: 10/08/2017

## simulate toy data
include("igp_mix.jl")
srand(09142017)
n = 60
n_iter = 1000
K = 3
c = rand(1:K, n)
update_ix = rand(1:n)
alpha = 1.0
thetas = Dict{Int64, KernelParam}()

a = GPHyper(
  Distributions.Logistic(-8, 1.5),
  Distributions.Logistic(0, 0.5),
  Distributions.Logistic(0, 0.5)
)

## consider real microbiome series
alpha = 0.15
y = readcsv("data/raw_data.csv")[:]
y += 0.01 * rand(length(y))
x = collect(linspace(0, 1, length(y)))[:, :]
MixGPSampler(x, y, alpha, a, "data/samples/", n_iter)

states = read_states(
  "data/samples/thetas.csv",
  "data/samples/c.csv"
)

x_new = collect(minimum(x):0.005:maximum(x))[:, :]
posteriors = mix_posteriors(x_new, states)
write_posteriors("data/posteriors.csv", x_new, posteriors)
writecsv("data/data.csv", [zeros(length(y)) x y])
