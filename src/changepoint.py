import os
import argparse
import numpy as np
import BASIC_changepoint as BASIC

parser = argparse.ArgumentParser(description="Apply the BASIC model.")
parser.add_argument(
    "--dir",
    type=str,
    default=os.path.join("..", "data", "changepoint"),
    help="Directory to read / write from"
)
parser.add_argument(
    "--csv",
    type=str,
    default="abt.csv",
    help="CSV name with the data to apply BASIC to"
)
parser.add_argument(
    "--iter",
    type=int,
    default=200,
    help="Number of iterations to run BASIC"
)
parser.add_argument(
    "--model",
    type=str,
    default="normal_var",
    help="BASIC likelihood model specifier"
)
args = parser.parse_args()

###############################################################################
#                               helper functions                              #
###############################################################################
def reshape_samples(samples):
    """
    Reshape Sample Dictionary

    To save results to a file, we convert the original samples dictionary into
    a "tall" array. The three columns are iteration, changepoint location, and
    rows associated with that location given that iteration.

    Arguments:
      samples [dictionary]: The output of a call to BASIC.MCMC_sample()
    """
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
X = np.loadtxt(
    os.path.join(args.dir, args.csv),
    delimiter=","
)

[samples, model_params, pi_q, q_vals] = BASIC.MCMC_sample(
    X,
    args.model,
    sample_iters=args.iter,
    MCEM_schedule=[10, 20, 40, 60, 100]
)

## writing outputs to file
np.savetxt(
    os.path.join(args.dir, "samples.csv"),
    reshape_samples(samples),
    delimiter=","
)

np.savetxt(
    os.path.join(args.dir, "params.csv"),
    np.vstack((
        model_params.keys(),
        model_params.values()
    )),
    fmt = "%s"
)

np.savetxt(os.path.join(args.dir, "q_vals.csv"), q_vals)
np.savetxt(os.path.join(args.dir, "pi_q.csv"), pi_q)
