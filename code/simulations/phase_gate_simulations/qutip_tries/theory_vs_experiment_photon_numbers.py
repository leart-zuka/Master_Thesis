import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Experimental data
photon_numbers_exp = [1.005e-5, 0.019, 0.0992, 0.2005, 0.4969, 0.9993, 4.3633]
photon_numbers_exp_err = [5.403e-8, 0.0001, 0.0002, 0.0002, 0.0007, 0.0012, 0.0236]
overlap_exp = [48.48, 85.03, 88.15, 85.30, 81.68, 80.42, 63.90]
overlap_exp_err = [11.87, 2.52, 0.86, 0.66, 0.52, 0.42, 0.20]

# Convert percentages to fractions
overlap_exp = np.array(overlap_exp) / 100
overlap_exp_err = np.array(overlap_exp_err) / 100

# Theoretical data from CSV
df = pd.read_csv("sim_data_photon_number_scan_test.csv")
photon_number_th = df["Photon Numbers"].values
overlap_th = df["Overlaps"].values

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(photon_number_th, overlap_th, "-", label="Theory")

df = pd.read_csv("sim_data_photon_number_scan_photon_dim_4.csv")
photon_number_th = df["Photon Numbers"].values
overlap_th = df["Overlaps"].values

ax.plot(photon_number_th, overlap_th, "-", label="Theory photon dim 4")

ax.errorbar(
    photon_numbers_exp,
    overlap_exp,
    xerr=photon_numbers_exp_err,
    yerr=overlap_exp_err,
    fmt="o",
    color="C1",
    capsize=3,
    label="Experiment",
)

ax.set_xlabel(r"Mean photon number $\bar{n}$")
ax.set_ylabel("Population Overlap")
ax.set_xscale("log")
ax.legend()
fig.tight_layout()
# plt.savefig("./plots/overlap_comparisson.pdf")
# plt.savefig("./plots/overlap_comparisson.svg")
plt.show()
