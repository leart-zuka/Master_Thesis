import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# ── Style setup (matching thesis figures) ────────────────────────────
mpl.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)

# ── Load data ────────────────────────────────────────────────────────
df = pd.read_csv("sim_data_photon_number_scan.csv", index_col=0)
n_bar = df["Photon Numbers"].values
overlap = df["Overlaps"].values

# Drop the n_bar = 0 point for log-scale plotting
mask = n_bar > 0
n_bar_plot = n_bar[mask]
overlap_plot = overlap[mask]

# ── Find optimum ─────────────────────────────────────────────────────
idx_max = np.argmax(overlap_plot)
n_opt = n_bar_plot[idx_max]
f_opt = overlap_plot[idx_max]

# ── Plot ─────────────────────────────────────────────────────────────
color_main = "#2060B0"

fig, ax = plt.subplots(figsize=(7, 4.5))

ax.semilogx(n_bar_plot, overlap_plot, color=color_main, lw=1.6)

# Vertical line at optimum
ax.axvline(x=n_opt, color="0.4", lw=0.9, ls="--", zorder=0)

# Reference lines
ax.axhline(y=0.5, color="0.75", lw=0.6, ls="--", zorder=0)

ax.set_xlabel(r"Mean photon number $\bar{n}$")
ax.set_ylabel("CNOT overlap")
ax.set_xlim(n_bar_plot[0], n_bar_plot[-1])
ax.set_ylim(0.45, 0.90)

fig.tight_layout()
fig.savefig("photon_number_scan.svg", bbox_inches="tight")
fig.savefig("photon_number_scan.pdf", bbox_inches="tight")
fig.savefig("photon_number_scan.png", dpi=200, bbox_inches="tight")
print(f"Optimum: n_bar = {n_opt:.2f}, Overlap = {f_opt:.4f}")
print("Done.")
