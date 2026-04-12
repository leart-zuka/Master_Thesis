import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Experimental data
photon_numbers_exp = [0.019, 0.0992, 0.2005, 0.4969, 0.9993, 4.3633]
photon_numbers_exp_err = [0.0001, 0.0002, 0.0002, 0.0007, 0.0012, 0.0236]
overlap_exp = [85.03, 88.15, 85.30, 78.11, 80.75, 65.49]
overlap_exp_err = [2.52, 0.86, 0.66, 0.65, 0.39, 0.36]

# Convert percentages to fractions
overlap_exp = np.array(overlap_exp) / 100
overlap_exp_err = np.array(overlap_exp_err) / 100

# Theoretical data from CSV
df = pd.read_csv("sim_data.csv")
photon_numbers_th = df["Photon Numbers"].values
fidelities_th = df["Overlaps"].values

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(photon_numbers_th, fidelities_th, "-", color="C0", label="Theory")
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
# ax.set_xscale("log")
# ax.set_xlim(1e-2, 2)
ax.set_ylim(0.5, 1.0)
ax.legend()
fig.tight_layout()
plt.savefig("./plots/overlap_comparisson.pdf")
plt.savefig("./plots/overlap_comparisson.svg")
plt.show()
