import glob
import numpy as np

OUTDIR = "out"


def get_combined_data(files):
    """Combine files into a single file"""
    data = []
    for f in files:
        data.append(np.loadtxt(f, skiprows=1))
    if len(data) == 0:
        return np.array([])
    return np.vstack(data)


def main():
    for d in [1, 20, 100]:
        data = get_combined_data(glob.glob(f"{OUTDIR}/*_d{d}_v*.dat"))
        if data.shape[0] == 0:
            continue
        print(f"Saving combined data for d={d} ({data.shape[0]} samples)")
        np.savetxt(f"{OUTDIR}/lartillot_d{d}_combined.dat", data, header="lnz lnz_err")


if __name__ == "__main__":
    main()
