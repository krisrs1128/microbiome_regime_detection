import numpy as np
import BASIC_changepoint as BASIC

###############################################################################
#                               helper functions                              #
###############################################################################
def reshape_samples(samples):
    output_matrix = np.empty((0, 3))
    for iter in range(len(samples)):
        for changept in samples[iter].keys():
            cur_rows = samples[iter][changept]
            n_cur = len(cur_rows)

            new_entry = np.hstack((
                iter * np.ones((n_cur, 1)),
                changept * np.ones((n_cur, 1)),
                np.matrix(cur_rows).reshape(n_cur, 1)
            ))
            output_matrix = np.vstack((
                output_matrix,
                new_entry
            ))
    return output_matrix

###############################################################################
#                               analysis script                                #
###############################################################################
## Gaussian mean / variance model
X = np.loadtxt("../data/changepoint/abt_asinh.csv", delimiter=",")

[samples, model_params, pi_q, q_vals] = BASIC.MCMC_sample(
    X,
    'normal_var',
    sample_iters=200,
    MCEM_schedule=[10, 20, 40, 60, 100]
)

np.savetxt(
    "../data/changepoint/samples.csv",
    reshape_samples(samples),
    delimiter=","
)
