import os
import numpy as np
import bilby
import logging
from bilby.core.likelihood import GaussianLikelihood
import pandas as pd

np.random.seed(42)
logging.getLogger('bilby').setLevel(logging.CRITICAL)


def get_priors():
    return bilby.core.prior.PriorDict(dict(
        m=bilby.core.prior.Uniform(0, 5, "m"),
        c=bilby.core.prior.Uniform(-2, 2, "c")
    ))


def generate_data(true_vals, noise_sigma, N=40):
    time = np.linspace(0, 10, N)
    noise = np.random.normal(0, noise_sigma, N)
    y = model(time, **true_vals) + noise
    return time, y


def model(time, m, c):
    return time * m + c


class GaussianLikelihoodWithNoiseLikelihood(GaussianLikelihood):
    def noise_log_likelihood(self):
        """Normally the GaussianLikelihood class returns a nan for the noise log likelihood."""
        return np.sum(- (self.y / self.sigma) ** 2 / 2 - np.log(2 * np.pi * self.sigma ** 2) / 2)


def run_inference():
    true_vals = dict(m=0.5, c=0.2)
    noise_sigma = 1
    time, data = generate_data(true_vals, noise_sigma)

    likelihood = GaussianLikelihood(time, data, model, sigma=noise_sigma)
    likelihood_w_noise_lnl = GaussianLikelihoodWithNoiseLikelihood(time, data, model, sigma=noise_sigma)

    sampling_kwgs = dict(
        priors=get_priors(),
        sampler="dynesty",
        nlive=1000,
        injection_parameters=true_vals,
        outdir='out',
        clean=True,
        verbose=False,
        plot=True,
        save=False
    )

    print("STARTING INFERENCE")
    result_no_noise_lnl = bilby.run_sampler(**sampling_kwgs, likelihood=likelihood, label='no_noise_lnl')
    result_with_noise_lnl = bilby.run_sampler(**sampling_kwgs, likelihood=likelihood_w_noise_lnl, label='with_noise_lnl')

    # print a table of the results
    for r in [result_no_noise_lnl, result_with_noise_lnl]:
        print(f"Result for {r.label}")
        print(f"LnZ: {r.log_evidence:.2f}")
        print(f"noise LnZ: {r.log_noise_evidence:.2f}")
        print(f"lnBF: {r.log_bayes_factor:.2f}")
        print("\n")


if __name__ == "__main__":
    run_inference()