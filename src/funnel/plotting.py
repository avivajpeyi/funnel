import numpy as np
import pandas as pd
from .fi_core import get_fi_lnz_list

import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 20})


def plot_fi_evidence_results(
        posterior_samples: pd.DataFrame = pd.DataFrame, sampling_lnz: np.array = [], r_vals: np.array = np.array([]),
        num_ref_params: int = 10, plot_inividual_lnzs: bool = False,
        lnzs=np.array([]),
):
    if len(lnzs) == 0:
        lnzs, r_vals = get_fi_lnz_list(posterior_samples, r_vals, num_ref_params)
    lnz_quants = np.nanquantile(lnzs, [0.16, 0.5, 0.84], axis=0)

    # PLOT Nested Sampling LnZ
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    if len(sampling_lnz) > 0:
        ax.fill_between(
            x=r_vals, y1=min(sampling_lnz), y2=max(sampling_lnz),
            color='tab:blue', interpolate=True, alpha=.5, label="Sampling LnZ",
            zorder=10
        )
    ax.set_xlim(min(r_vals), max(r_vals))
    ax.set_xlabel(r"FI $R$")
    ax.set_ylabel(r"$\ln{\mathcal{Z}}$")
    ax.set_xscale('log')

    if plot_inividual_lnzs:
        for lnz in lnzs:
            ax.plot(r_vals, lnz, color="tab:green", alpha=.05, zorder=-1)
    else:
        ax.fill_between(
            x=r_vals, y1=lnz_quants[2], y2=lnz_quants[0],
            color='tab:green', interpolate=True, alpha=.25, zorder=0, lw=0, label=r"90% CI"
        )
        ax.plot(r_vals, lnz_quants[1], color="tab:green", alpha=1, zorder=1)
    ax.plot([], [], label="FI LnZ", color="tab:green", alpha=1, zorder=1)

    ax.legend(loc=(1.1, 0.5), frameon=False)
    plt.tight_layout()
    return fig
