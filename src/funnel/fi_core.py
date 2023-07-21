import numpy as np
import bilby
from tqdm.auto import tqdm
import matplotlib.pyplot as plt


def ln_posterior_approx_for_sample(posterior_samples, ref_parm, r):
    diff_from_ref = posterior_samples - ref_parm
    sin_diff = np.sin(r * diff_from_ref)

    sin_diff_divided_diff = sin_diff / diff_from_ref
    res_prod = np.prod(sin_diff_divided_diff, axis=1)
    sum_res = np.sum(res_prod)

    sum_res = np.sum(np.prod(sin_diff / diff_from_ref, axis=1))

    # log_terms = np.log(sin_diff / diff_from_ref)
    # max_val = log_terms.max(axis=1)
    # res_2 = max_val + np.log(np.sum(np.exp(log_terms - max_val[:, np.newaxis]), axis=1))

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
    r_vals = np.geomspace(1e-3, 1e5, 1000)

    sampling_lnz_up = result.log_evidence + result.log_evidence_err
    sampling_lnz_low = result.log_evidence - result.log_evidence_err

    # PLOT Nested Sampling LnZ
    fig, ax = plt.subplots(1, 1)

    ax.fill_between(
        x=r_vals, y1=sampling_lnz_up, y2=sampling_lnz_low,
        color='tab:blue', interpolate=True, alpha=.5, label="Nested Samp. LnZ"
    )
    ax.set_xlim(min(r_vals), max(r_vals))
    ax.set_xlabel(r"FI $R$")
    ax.set_ylabel(r"$\ln{\mathcal{Z}}$")
    ax.set_xscale('log')

    if len(ref_parms) == 0:
        ref_parms = [result.injection_parameters]

    for i, ref_parm in tqdm(enumerate(ref_parms)):
        ref_parm = {k: ref_parm[k] for k in result.search_parameter_keys}
        ref_lnpri, ref_lnl = get_refernce_lnprior_lnlikelihood(priors, likelihood, ref_parm)

        fi_kwargs = dict(
            posterior_samples=result.posterior[result.search_parameter_keys].values,
            ref_parm=np.array([*ref_parm.values()]),
            reference_lnprior=ref_lnpri,
            reference_lnlikelihood=ref_lnl,
        )
        ln_evidences = [fi_ln_evidence(**fi_kwargs, r=ri) for ri in r_vals]
        print(ln_evidences)
        if i == 0:
            ax.plot(r_vals, ln_evidences, label="FI LnZ", color="tab:green", alpha=0.70)
        else:
            ax.plot(r_vals, ln_evidences, color="tab:green", alpha=0.70)

    ax.legend(loc=(1.1, 0.8), frameon=False)
    plt.tight_layout()
    return fig
