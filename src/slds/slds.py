"""
Script to apply SLDS to a single series

author: sankaran.kris@gmail.com
date: 10/13/2017
"""

import numpy as np
import pyslds.models as slds

###############################################################################
## fitting functions
###############################################################################
def slds_fit(y, n_iter, K=3, D_latent=1, outdir="."):
    outpaths = [
        os.path.join(outdir, "emission.csv"),
        os.path.join(outdir, "dynamics.csv"),
        os.path.join(outdir, "stateseq.csv"),
        os.path.join(outdir, "rvs.csv")
    ]
    for o in outpaths:
        if os.path.exists(o):
            os.remove(o)

    inputs = np.zeros((y.shape[1], 0))
    for i in range(y.shape[0]):
        model = slds.DefaultSLDS(K, 1, D_latent, 0)
        model.add_data(y[i, ].reshape(y.shape[1], 1))

        for iter in range(n_iter):
            model.resample_model()
            if (iter % 5 == 0):
                print(
                    "\t".join([str(s) for s in ["sequence", i, "iter ", iter]])
                )

                write_parameters(i, iter, model.emission_distns, outpaths[0], ["C", "R"])
                write_parameters(i, iter, model.dynamics_distns, outpaths[1], ["A", "Q"])
                write_states(i, iter, model.stateseqs, outpaths[2])
                write_rvs(i, iter, model.rvs(), outpaths[3])


def write_states(i, iter, stateseqs, outpath):
    with open(outpath, "a") as f:
        cols = [i, iter] + stateseqs[0].tolist()
        f.write(",".join([str(x) for x in cols]) + "\n")


def write_rvs(i, iter, rvs, outpath):
    with open(outpath, "a") as f:
        list_rvs = sum(rvs.tolist(), [])
        for l, v in enumerate(list_rvs):
            cols = [i, iter, l, v]
            f.write(",".join([str(x) for x in cols]) + "\n")


def write_parameters(i, iter, reglist, outpath, pnames):
    for k, reg in enumerate(reglist):
        for j, theta in enumerate(reg.parameters):
            list_params = sum(theta.tolist(), [])
            with open(outpath, "a") as f:
                for l, v in enumerate(list_params):
                    cols = [i, iter, k, pnames[j], l, v]
                    f.write(",".join([str(x) for x in cols]) + "\n")


###############################################################################
##  load data and run models
###############################################################################
y = np.loadtxt("../../data/slds/abt.csv", delimiter=",")
slds_fit(y[1:10, :], 200, 3, 1, "../../data/slds/")
