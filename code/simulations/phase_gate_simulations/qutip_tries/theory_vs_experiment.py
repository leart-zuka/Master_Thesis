import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Experimental data
photon_numbers_exp = [0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.6]
overlap_exp = [78.05, 86.75, 91.33, 90.85, 86.39, 86.88, 82.52, 84.84]
overlap_exp_err = [4.58, 2.11, 1.86, 1.25, 1.05, 2.40, 1.51, 0.84]

# Convert percentages to fractions
overlap_exp = np.array(overlap_exp) / 100
overlap_exp_err = np.array(overlap_exp_err) / 100

# Theoretical data from CSV
df = pd.read_csv("sim_data.csv")
photon_numbers_th = df["Photon Numbers"].values
fidelities_th = df["Fidelities"].values

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(photon_numbers_th, fidelities_th, "-", color="C0", label="Theory")
ax.errorbar(
    photon_numbers_exp,
    overlap_exp,
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
plt.savefig("fidelity_comparison.pdf")
plt.show()
