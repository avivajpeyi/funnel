import numpy as np
import pandas as pd
from .utils import get_post_mask
from .logger import logger
from tqdm.auto import trange
import os


def fi_ln_posterior(posterior_samples: np.ndarray, reference_sample: np.array, r: float):
    """
    Returns the approx posterior probability at the reference parameter value.
    The approximation is based on the 'density estimation' method described in
    [Rotiroti et al., 2018](https://link.springer.com/article/10.1007/s11222-022-10131-0).

    :param posterior_samples:np.ndarray: Array of posterior samples [n_samples, n_dim]
    :param reference_sample:np.array: A reference parameter value [n_dim] (Not present in the posterior)
    :param r:float: A scaling factor
    :return: The log of the approximated posterior
    """

    diff_from_ref = posterior_samples - reference_sample

    sin_diff = np.sin(r * diff_from_ref)
    sum_res = np.nansum(np.nanprod(sin_diff / diff_from_ref, axis=1))
    n_samp, n_dim = posterior_samples.shape
    const = 1 / (n_samp * np.power(np.pi, n_dim))
    if sum_res < 0:
        return np.nan
    return np.log(sum_res * const)


def fi_ln_evidence(posterior_samples: np.ndarray, ref_samp: np.array, r: float, ref_lnpri: float, ref_lnl: float):
    """
    Returns the approx log-evidence of some posterior samples (using a reference parameter value).
    The approximation is based on the 'density estimation' method described in
    [Rotiroti et al., 2018](https://link.springer.com/article/10.1007/s11222-022-10131-0).

    :param posterior_samples:np.ndarray: Array of posterior samples [n_samples, n_dim]
    :param ref_samp:np.array: A reference parameter value [n_dim] (Not present in the posterior)
    :param r:float: A scaling factor
    :param ref_lnpri:float: The log of the reference prior
    :param ref_lnl:float: The log of the reference likelihood
    :return: The log of the approximated log-evidence
    """
    approx_ln_post = fi_ln_posterior(posterior_samples, ref_samp, r)
    return ref_lnpri + ref_lnl - approx_ln_post


def get_fi_lnz_list(
        posterior_samples: pd.DataFrame, r_vals: np.array = [], num_ref_params: int = 10,
        cache_fn=""
):
    if os.path.exists(cache_fn):
        data = np.load(cache_fn)
        return data['lnzs'], data['r_vals']

    if len(r_vals) == 0:
        r_vals = np.geomspace(1e-3, 1e10, 2000)

    if num_ref_params > len(posterior_samples):
        num_ref_params = len(posterior_samples)

    # unpacking posterior data
    ln_pri = posterior_samples['log_prior'].values
    ln_lnl = posterior_samples['log_likelihood'].values
    post = posterior_samples[posterior_samples.columns.drop(["log_prior", "log_likelihood"])].values

    logger.info(f"Calculating FI LnZ with {num_ref_params} reference points "
                f"and a posterior of size: {post.shape}")

    # randomly select reference points
    ref_idx = np.random.choice(len(post), num_ref_params, replace=False)

    lnzs = np.zeros((num_ref_params, len(r_vals)))
    median_lnzs = np.zeros(num_ref_params)
    med_ = 0

    with trange(num_ref_params, desc="FI LnZ", postfix=f"FI LnZ: {med_}") as pbar:
        for i in pbar:
            refi = ref_idx[i]
            med_ = np.nanmedian(median_lnzs[:i]) if i > 0 else 0

            post_mask = get_post_mask(post, refi)
            fi_kwargs = dict(
                posterior_samples=post[post_mask],
                ref_samp=post[refi],
                ref_lnpri=ln_pri[refi],
                ref_lnl=ln_lnl[refi],
            )
            lnzs[i] = np.array([fi_ln_evidence(**fi_kwargs, r=ri) for ri in r_vals])
            median_lnzs[i] = np.nanmedian(lnzs[i])
            pbar.set_postfix_str(f"FI LnZ: {med_:.2f}")
            pbar.update()

    if cache_fn:
        np.savez(cache_fn, lnzs=lnzs, r_vals=r_vals)

    return lnzs, r_vals
