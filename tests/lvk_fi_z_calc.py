import h5py
import numpy as np
import warnings
import logging
import pandas as pd
from collections import namedtuple

from funnel.plotting import plot_fi_evidence_results

logging.getLogger('bilby').setLevel(logging.CRITICAL)

# re-defining plotting defaults
from matplotlib import rcParams



# load the LVK posterior
FPATH = '../../NR_DATA/IGWN-GWTC2p1-v2-GW150914_095045_PEDataRelease_mixed_cosmo.h5'

LVK_data = namedtuple("LVK_data", "posterior, lnz, lnz_err")

def load_lvk_data(fpath):
    with h5py.File(fpath, 'r') as f:
        sampling_params = list(f['C01:IMRPhenomXPHM/priors/analytic'].keys())
        print(f['C01:IMRPhenomXPHM/priors/calibration'].keys())
        sampler_data = f['C01:IMRPhenomXPHM/meta_data/sampler']
        lnz, lnz_err =  sampler_data['ln_evidence'][0], sampler_data['ln_evidence_error'][0]
        post = f['C01:IMRPhenomXPHM']['posterior_samples'][()]
        post = pd.DataFrame({name: post[name][:] for name in post.dtype.names})
        post = post[sampling_params + ['log_likelihood', 'log_prior']]
        # post drop any columns with one unique value
        post = post.loc[:, post.nunique() > 1]
    return LVK_data(post, lnz, lnz_err)

GW150914_data = load_lvk_data(FPATH)

for c in GW150914_data.posterior.columns:
    print(c)


post = GW150914_data.posterior.sample(10000)
fig = plot_fi_evidence_results(
    post,
    # sampling_lnz=[GW150914_data.lnz+GW150914_data.lnz_err, GW150914_data.lnz-GW150914_data.lnz_err],
    num_ref_params=10,
    r_vals=np.geomspace(50, 2000, 10),
)
fig.show()
