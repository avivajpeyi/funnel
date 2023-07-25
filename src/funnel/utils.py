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
def get_post_mask(post: np.ndarray, idx: int, tol=1e-5) -> np.array:
    """Get the posterior mask after excluding the idx-th row (and its close neighbors)"""
    ref_samp = post[idx]
    distances = np.linalg.norm(post - ref_samp, axis=1)
    mask = distances > tol
    num_close = sum(~mask)

    if num_close > 3:
        logger.warning(f"Found {num_close} posterior samples that are close to ref_samp (tol={tol}).")

    return mask
