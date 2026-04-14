import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Experimental data
transmission_exp = [0.7738, 0.3182, 0.2247]
transmission_exp_err = [0.0013, 0.0010, 0.0004]
overlap_exp = [79.02, 84.78, 88.15]
overlap_exp_err = [1.05, 0.92, 0.86]

# Convert percentages to fractions
overlap_exp = np.array(overlap_exp) / 100
overlap_exp_err = np.array(overlap_exp_err) / 100

# Theoretical data from CSV
df = pd.read_csv("sim_data_bw_scan_scan.csv")
transmissions_th = df["Transmissions"].values
overlaps_th = df["Overlaps"].values

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(transmissions_th, overlaps_th, "-", color="C0", label="Theory")
ax.errorbar(
    transmission_exp,
    overlap_exp,
    xerr=transmission_exp_err,
    yerr=overlap_exp_err,
    fmt="o",
    color="C1",
    capsize=3,
    label="Experiment",
)

ax.set_xlabel(r"Transmission for the V polarization")
ax.set_ylabel("Population Overlap")
ax.legend()
fig.tight_layout()
plt.savefig("./plots/overlap_comparisson_bw.pdf")
plt.savefig("./plots/overlap_comparisson_bw.svg")
plt.show()
