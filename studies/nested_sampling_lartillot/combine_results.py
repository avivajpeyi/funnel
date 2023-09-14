import glob

import numpy as np

OUTDIR = "out"


def get_combined_data(files):
    """Combine files into a single file"""
    data = []
    for f in files:
        # if file just has 1 line, skip it
        if sum(1 for line in open(f)) == 1:
            continue
        data.append(np.loadtxt(f, skiprows=1))
    if len(data) == 0:
        return np.array([])
    return np.vstack(data)


def main():
    data = []
    for d in [1, 20, 100]:
        data_di = get_combined_data(glob.glob(f"{OUTDIR}/*_d{d}_v*.dat"))
        if data_di.shape[0] == 0:
            continue
        data_di = np.hstack([np.ones((data_di.shape[0], 1)) * d, data_di])
        data.append(data_di)
    data = np.vstack(data)
    print(f"Saving combined data for d={d} ({data.shape[0]} samples)")
    np.savetxt(
        f"{OUTDIR}/nested_sampling_lnzs.dat",
        data,
        header="dim ns_lnz ns_lnz_err",
        fmt=["%03d", "%.5e", "%.5e"],
    )


if __name__ == "__main__":
    main()
