import os

import pytest
import numpy as np
import bilby
import matplotlib.pyplot as plt
import shutil

from funnel.fi_core import get_refernce_lnprior_lnlikelihood, plot_fi_evidence_results, fi_ln_evidence

CLEAN = False

np.random.seed(42)

@pytest.fixture(scope="module")
def bilby_objects():

    def model(time, m, c):
        return time * m + c

    outdir, label = 'out', 'line'
    if os.path.exists(outdir) and CLEAN:
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Now we define the injection parameters which we make simulated data with
    injection_parameters = dict(m=0.5, c=0.2)

    sampling_frequency = 10
    time_duration = 10
    time = np.arange(0, time_duration, 1 / sampling_frequency)
    N = len(time)
    sigma = np.random.normal(1, 0.01, N)
    data = model(time, **injection_parameters) + np.random.normal(0, sigma, N)

    likelihood = bilby.likelihood.GaussianLikelihood(time, data, model, sigma)

    priors = bilby.core.prior.PriorDict(dict(
        m=bilby.core.prior.Uniform(0, 5, "m"),
        c=bilby.core.prior.Uniform(-2, 2, "c")
    ))

    # We quickly plot the data to check it looks sensible
    if not os.path.exists(outdir):
        fig, ax = plt.subplots()
        ax.plot(time, data, "o", label="data")
        ax.plot(time, model(time, **injection_parameters), "--r", label="signal")
        ax.set_xlabel("time")
        ax.set_ylabel("y")
        ax.legend()
        plt.savefig(f"{outdir}/data.png")

    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=priors,
        sampler="dynesty",
        nlive=1500,
        outdir=outdir,
        label=label,
        injection_parameters=injection_parameters,
    )
    fig = result.plot_corner()
    fig.suptitle(f"Nested Sampling LnZ = {result.log_evidence:.2f}+/-{result.log_evidence_err:.2f}", y=1.1)
    fig.savefig(f"{outdir}/corner.png")

    return result, priors, likelihood


def test_fi_integration_one_val(bilby_objects):
    result, priors, likelihood = bilby_objects

    ref_lnpri, ref_lnl = get_refernce_lnprior_lnlikelihood(priors, likelihood, result.injection_parameters)

    fi_lnz_test = fi_ln_evidence(
        posterior_samples=result.posterior[['m', 'c']].values,
        ref_parm=np.array([*result.injection_parameters.values()]),
        r=100,
        reference_lnprior=ref_lnpri,
        reference_lnlikelihood=ref_lnl,
    )
    print(f"Test FI LnZ (R=100) = {fi_lnz_test:.2f}")
    assert np.isfinite(fi_lnz_test)


def test_fi_integration_plot(bilby_objects):
    result, priors, likelihood = bilby_objects
    # fig = plot_fi_evidence_results(result, priors, likelihood)
    # fig.suptitle("FI Evidence @ True Injection", y=1.1)
    # fig.tight_layout()
    # fig.savefig(f"{result.outdir}/fi_evidence_at_true_inj.png")
    n_samp = 100
    fig = plot_fi_evidence_results(result, priors, likelihood, result.posterior.sample(n_samp).to_dict('records'))
    fig.suptitle(f"FI Evidence @ {n_samp} Posterior Samples", y= 1.1)
    fig.tight_layout()
    fig.savefig(f"{result.outdir}/fi_evidence_{n_samp}_posterior_samples_.png")