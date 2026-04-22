import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)

# ── EXPERIMENTAL DATA ──
bw_transmission = [0.7738, 0.3182, 0.2247]
bw_transmission_err = [0.0013, 0.0010, 0.0004]
bw_overlap = np.array([79.02, 84.78, 88.15]) / 100
bw_overlap_err = np.array([1.05, 0.92, 0.86]) / 100

pn_numbers = [1.005e-5, 0.019, 0.0992, 0.2005, 0.4969, 0.9993, 4.3633]
pn_numbers_err = [5.403e-8, 0.0001, 0.0002, 0.0002, 0.0007, 0.0012, 0.0236]
pn_overlap = np.array([48.48, 85.03, 88.15, 85.30, 81.68, 80.42, 63.90]) / 100
pn_overlap_err = np.array([11.87, 2.52, 0.86, 0.66, 0.52, 0.42, 0.20]) / 100

# ── THEORY DATA ──
df_bw = pd.read_csv("sim_data_bw_scan_scan.csv")
bw_trans_th = df_bw["Transmissions"].values
bw_overlap_th = df_bw["Overlaps"].values

df_pn = pd.read_csv("sim_data_photon_number_scan.csv")
pn_number_th = df_pn["Photon Numbers"].values
pn_overlap_th = df_pn["Overlaps"].values

# ── COLORS ──
theory_color = "#2b6ca3"
exp_color = "#c97a2a"

# ── FIGURE ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# ── LEFT: Brewster window scan ──
ax1.plot(
    bw_trans_th,
    bw_overlap_th,
    "-",
    color=theory_color,
    linewidth=1.5,
    label="Simulation",
    zorder=2,
)
ax1.errorbar(
    bw_transmission,
    bw_overlap,
    xerr=bw_transmission_err,
    yerr=bw_overlap_err,
    fmt="o",
    color=exp_color,
    markersize=6,
    capsize=3,
    capthick=1,
    elinewidth=1,
    label="Experiment",
    zorder=3,
)
ax1.set_xlabel(r"$V$-mode transmission")
ax1.set_ylabel("Population overlap")
ax1.set_xlim(0.0, 1.0)
ax1.set_ylim(0.45, 1.0)
ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
ax1.legend(frameon=True, fancybox=False, edgecolor="#cccccc")

# ── RIGHT: Photon number scan ──
ax2.plot(
    pn_number_th,
    pn_overlap_th,
    "-",
    color=theory_color,
    linewidth=1.5,
    label="Simulation",
    zorder=2,
)
ax2.errorbar(
    pn_numbers,
    pn_overlap,
    xerr=pn_numbers_err,
    yerr=pn_overlap_err,
    fmt="o",
    color=exp_color,
    markersize=6,
    capsize=3,
    capthick=1,
    elinewidth=1,
    label="Experiment",
    zorder=3,
)
ax2.set_xlabel(r"Mean photon number $\bar{n}$")
ax2.set_ylabel("Population overlap")
ax2.set_xscale("log")
ax2.set_ylim(0.35, 1.0)
ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
ax2.legend(frameon=True, fancybox=False, edgecolor="#cccccc")

# ── PANEL LABELS ──
ax1.text(
    -0.15,
    1.05,
    "(a)",
    transform=ax1.transAxes,
    fontsize=12,
    fontweight="bold",
    va="top",
)
ax2.text(
    -0.15,
    1.05,
    "(b)",
    transform=ax2.transAxes,
    fontsize=12,
    fontweight="bold",
    va="top",
)

fig.tight_layout(w_pad=3)
# plt.show()
plt.savefig("sim_comparison_plots.svg", bbox_inches="tight")
plt.savefig("sim_comparison_plots.png", bbox_inches="tight", dpi=250)
print("Done")
