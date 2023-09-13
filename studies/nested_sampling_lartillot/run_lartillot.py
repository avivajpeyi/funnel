import bilby
from funnel.fi_core import get_fi_lnz_list
import numpy as np
import pandas as pd
import argparse
from tqdm.auto import tqdm, trange
from contextlib import redirect_stdout
import io

import logging

logging.getLogger("bilby").setLevel(logging.ERROR)


class LartillotLikelihood(bilby.Likelihood):
    def __init__(self, dim, v):
        self.dim = dim
        self.v = v
        super().__init__(parameters={f"x{i}": None for i in range(dim)})

    def log_likelihood(self):
        x = np.array(list(self.parameters.values()))
        lnl = -np.power(x, 2) / (2 * self.v)
        if self.dim == 1:
            return lnl[0]
        else:
            return np.sum(lnl)


def get_prior(n_dim):
    pri = bilby.core.prior.PriorDict()
    for i in range(n_dim):
        pri[f"x{i}"] = bilby.core.prior.TruncatedNormal(
            mu=0, sigma=1, minimum=-1000, maximum=1000
        )
    return pri


def true_lnz(v, dim):
    return (dim / 2) * (np.log(v) - np.log(1 + v))


def nested_sampling_lnz(v, dim):
    likelihood = LartillotLikelihood(dim, v)
    priors = get_prior(dim)
    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=priors,
        sampler="dynesty",
        nlive=1000,
        label=f"lartillot_dynesty_d{dim}_v{v}",
        clean=False,
        sample="rwalk",
    )
    return (result.log_evidence, result.log_evidence_err)


def parallel_tempering_lnz(v, dim):
    likelihood = LartillotLikelihood(dim, v)
    priors = get_prior(dim)
    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=priors,
        sampler="bilby_mcmc",
        ntemps=10,
        label=f"lartillot_ptmc_d{dim}_v{v}",
        clean=False,
        proposal_cycle="default",
        printdt=10,
        nsamples=2000,
    )
    return (result.log_evidence, result.log_evidence_err)


def __simulate_posterior(v, dim, nsamp=1000):
    likelihood = LartillotLikelihood(dim, v)
    priors = get_prior(dim)
    post = []
    scale = np.sqrt(v / (v + 1))
    for i in range(nsamp):
        params = {f"x{i}": np.random.normal(scale=scale) for i in range(dim)}
        likelihood.parameters.update(params)
        post.append(
            dict(
                log_likelihood=likelihood.log_likelihood(),
                log_prior=priors.ln_prob(likelihood.parameters),
                **params,
            )
        )
    return pd.DataFrame(post)


def fi_lnz(v, dim, nsamp=1000):
    lnzs, r_vals, samp = get_fi_lnz_list(
        __simulate_posterior(v, dim, nsamp=nsamp),
        r_vals=np.linspace(0.1, 100, 100),
        num_ref_params=1,
        weight_samples_by_lnl=True,
    )
    # only keep last 90% of lnzs
    lnzs = lnzs[:, -int(0.9 * nsamp) :]
    return np.median(lnzs), np.std(lnzs)


# make a main function that runs the nested_sampling_lnz given a v and dim from the command line
# save the LnZ to a file
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--v", type=float, default=0.01)
    parser.add_argument("-d", "--dim", type=int, default=1)
    parser.add_argument("-r", "--nrep", type=int, default=50)
    parser.add_argument("-s", "--seed", type=int, default=-1)
    args = parser.parse_args()
    if args.seed < 0:
        args.seed = np.random.randint(0, 1e5)
        np.random.seed(args.seed)

    print("Lartillot LnZ (v={args.v}, dim={args.dim}, seed={args.seed})")
    print(f"True LnZ: {true_lnz(args.v, args.dim):.2f}")
    outfile = f"lartillot_d{args.dim}_v{args.v}_seed{args.seed}.dat"
    pbar = tqdm(total=args.nrep)
    for _ in pbar:
        f = io.StringIO()
        with redirect_stdout(f):  # suppress output
            lnz, lnz_err = nested_sampling_lnz(args.v, args.dim)
        with open(outfile, "a") as f:
            f.write(f"{lnz} {lnz_err}\n")
        pbar.set_postfix_str(f"LnZ: {lnz:.2f} +/- {lnz_err:.2f}")


if __name__ == "__main__":
    main()
