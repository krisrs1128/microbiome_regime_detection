"""
Script to apply SLDS to a single series

author: sankaran.kris@gmail.com
date: 10/13/2017
"""

import numpy as np
import pyslds.models as slds

K = 3               # Number of discrete latent states
D_obs = 1           # Observed data dimension
D_latent = 2        # Latent state dimension
D_input = 0         # Exogenous input dimension

y = np.loadtxt("../../data/slds/abt.csv")
y = y.reshape(len(y), 1)
inputs = np.zeros((y.shape[0], 0))
model = slds.DefaultSLDS(K, D_obs, D_latent, D_input)
model.add_data(y)

# Run the Gibbs sampler
N_samples = 1000
def update(model):
    model.resample_model()
    print(model.stateseqs)
    return model.log_likelihood()

lls = [update(model) for _ in range(N_samples)]

## some other info we should probably track
model.emission_distns[0].parameters
model.dynamics_distns[0].parameters
model.trans_distn.trans_matrix


## could simulate from fitted models using code from
## https://github.com/mattjj/pylds/blob/12fe0d4a4b2dc30926e8b263ea1a4de034f65186/pylds/models.py
