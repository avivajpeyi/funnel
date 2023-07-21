import numpy as np
import bilby
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")


def ln_posterior_approx_for_sample(posterior_samples, ref_parm, r):
    diff_from_ref = posterior_samples - ref_parm
    sin_diff = np.sin(r * diff_from_ref)
    sum_res = np.sum(np.prod(sin_diff / diff_from_ref, axis=1))
    n_samp, n_dim = posterior_samples.shape
    const = 1 / (n_samp * np.power(np.pi, n_dim))
    return np.log(sum_res * const)


def fi_ln_evidence(posterior_samples, ref_parm, r, reference_lnprior, reference_lnlikelihood):
    approx_ln_post = ln_posterior_approx_for_sample(posterior_samples, ref_parm, r)
    return reference_lnprior + reference_lnlikelihood - approx_ln_post


def get_refernce_lnprior_lnlikelihood(prior, likelihood, ref_parm):
    likelihood.parameters.update(ref_parm)
    lnlike = likelihood.log_likelihood()
    lnprior = prior.ln_prob(ref_parm)
    return lnprior, lnlike


def plot_fi_evidence_results(result, priors, likelihood, ref_parms=[]):
    r_vals = np.geomspace(1e-3, 1e10, 2000)

    sampling_lnz_up = result.log_evidence + result.log_evidence_err
    sampling_lnz_low = result.log_evidence - result.log_evidence_err

    # PLOT Nested Sampling LnZ
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    ax.fill_between(
        x=r_vals, y1=sampling_lnz_up, y2=sampling_lnz_low,
        color='tab:blue', interpolate=True, alpha=.5, label="Nested Samp. LnZ",
        zorder=10
    )
    ax.set_xlim(min(r_vals), max(r_vals))
    ax.set_xlabel(r"FI $R$")
    ax.set_ylabel(r"$\ln{\mathcal{Z}}$")
    ax.set_xscale('log')

    if len(ref_parms) == 0:
        ref_parms = [result.injection_parameters]

    lnzs = []
    for i, ref_parm in tqdm(enumerate(ref_parms)):
        ref_parm = {k: ref_parm[k] for k in result.search_parameter_keys}
        ref_lnpri, ref_lnl = get_refernce_lnprior_lnlikelihood(priors, likelihood, ref_parm)

        post = result.posterior[result.search_parameter_keys].values
        ref_parm = np.array([*ref_parm.values()])

        # check if ref_parm is in posterior samples
        point_idx = np.argwhere(post == ref_parm)
        if len(point_idx) > 0:
            # remove ref_parm from posterior samples
            post = np.delete(post, point_idx[0][0], axis=0)

        fi_kwargs = dict(
            posterior_samples=post,
            ref_parm=ref_parm,
            reference_lnprior=ref_lnpri,
            reference_lnlikelihood=ref_lnl,
        )
        ln_evidences = [fi_ln_evidence(**fi_kwargs, r=ri) for ri in r_vals]
        lnzs.append(ln_evidences)

    lnzs = np.array(lnzs)
    lnz_quantiles = np.nanquantile(lnzs, [0.16, 0.5, 0.84], axis=0)
    lnz_up = lnz_quantiles[2]
    lnz_low = lnz_quantiles[0]
    lnz_med = lnz_quantiles[1]

    ax.plot(r_vals, lnz_med, label="FI LnZ", color="tab:green", alpha=1, zorder=1)
    ax.fill_between(
        x=r_vals, y1=lnz_up, y2=lnz_low,
        color='tab:green', interpolate=True, alpha=.25, zorder=0
    )
    for lnz in lnzs:
        ax.plot(r_vals, lnz, color="tab:green", alpha=.1, zorder=-1)

    ax.legend(loc=(1.1, 0.8), frameon=False)
    plt.tight_layout()
    return fig
