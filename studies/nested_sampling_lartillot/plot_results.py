import numpy as np
import matplotlib.pyplot as plt


# set fontsize to 20 for everything in rcPara
plt.rcParams.update({"font.size": 20})

nlive50_dat = np.loadtxt("lartillot_d1_v0.01_nlive50.dat", skiprows=1)
nlive1000_dat = np.loadtxt("lartillot_d1_v0.01_nlive1000.dat", skiprows=1)


bins = np.linspace(-2.6, -2, 10)

fig, ax = plt.subplots(2, 1, figsize=(5, 8))
ax[0].hist(nlive50_dat[:, 0], bins=bins, alpha=0.5, label="nlive=50", density=True)
ax[0].hist(nlive1000_dat[:, 0], bins=bins, alpha=0.5, label="nlive=1000", density=True)
ax[0].legend()
ax[0].set_xlabel("logZ")
ax[0].set_ylabel("p(logZ)")

# plot runtimes on y, label on x
runtimes = dict(
    nlive50=((2 * 60 + 20) / 50, np.std((2 * 60 + 20) / 50) / np.sqrt(50)),
    nlive1000=((9 * 60 + 46) / 23, np.std((9 * 60 + 46) / 23) / np.sqrt(23)),
)

ax[1].errorbar(
    [0, 1],
    [runtimes["nlive50"][0], runtimes["nlive1000"][0]],
    yerr=[1, 6],
    fmt="o",
)
ax[1].set_ylabel("runtime [s]")
ax[1].set_xticks([0, 1], ["nlive=50", "nlive=1000"])


plt.tight_layout()
plt.show()
