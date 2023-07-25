from bilby.gw.result import CBCResult
from funnel.plotting import plot_fi_evidence_results
import numpy as np

res = CBCResult.from_hdf5("injection_merged_result.hdf5")
lnz, lnz_err = res.log_evidence, res.log_evidence_err
post = res.posterior
parms = res.search_parameter_keys
post = post[parms + ['log_likelihood', 'log_prior']]


fig = plot_fi_evidence_results(
    post.sample(3000),
    # sampling_lnz=[lnz+lnz_err, lnz-lnz_err],
    num_ref_params=10,
    r_vals=np.geomspace(50, 10000, 100),
)
fig.show()