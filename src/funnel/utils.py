import pandas as pd
import numpy as np
from typing import Tuple
import time
from .logger import logger

LNPRI_LNL_COL_NAMES = ['log_prior', 'log_likelihood']


# make a decorator to time functions
def timefunct(func):
    def timed(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        runtime = end - start
        if runtime > 0.01:
            logger.debug(f"Runtime {func.__name__}: {end - start:.2f}s")
        return result

    return timed


@timefunct
def extract_ref_samp_from_post(post: pd.DataFrame) -> Tuple[np.ndarray, np.array, float, float]:
    """Extract ref param from posterior and return formatted data
    :return: post, ref_samp, ref_lnl, ref_lnpri
    """
    ln_pri_lnl = post[LNPRI_LNL_COL_NAMES].values
    post = post[post.columns.drop(LNPRI_LNL_COL_NAMES)].values
    idx = np.random.randint(0, len(post))
    ref_samp, ref_lnl, ref_pri = post[idx], ln_pri_lnl[idx, 0], ln_pri_lnl[idx, 1]

    # post without ref_samp_idx
    deletion_idx = get_idx_of_rows_with_similar_values(post, ref_samp, tol=1e-5)
    # mask out deletion_idx
    if len(deletion_idx) > 3:
        logger.warning(f"Found {len(deletion_idx)} posterior samples that are close to ref_samp (tol=1e-5).")
    post = post[~np.isin(np.arange(len(post)), deletion_idx)]
    return post, ref_samp, ref_lnl, ref_pri


def get_idx_of_rows_with_similar_values(arr: np.ndarray, sample_row: np.array, tol: float = 1e-10) -> np.ndarray:
    """Get the row idx in arr that are close to sample_row"""
    mask = np.isclose(arr, sample_row, rtol=tol)
    return np.where(mask.all(axis=1))[0]



def __check_ref_not_in_posterior(post, ref, tol=1e-5):
    close_rows = get_idx_of_rows_with_similar_values(post, ref, tol=tol)
    if len(close_rows) > 0:
        raise ValueError(f"Found {len(close_rows)} posterior samples that are close to ref_samp (tol={tol}).")
