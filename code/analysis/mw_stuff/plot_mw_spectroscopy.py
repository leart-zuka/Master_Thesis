import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("03_04_26_MWSpec_F1mf0_1.csv")

# Peak centers and band half-width (kHz)
peak_centers = [-100, 0, 100]
band_half_width = 15

# Yellow, green, light blue — for the three F=1,mF=0 transitions
colors = ["#f1c40f", "#2ecc71", "#00bcd4"]

fig, ax = plt.subplots(figsize=(11, 5.5))

# Colored vertical bands
for center, color in zip(peak_centers, colors):
    ax.axvspan(
        center - band_half_width,
        center + band_half_width,
        alpha=0.3,
        color=color,
        zorder=0,
    )

# Data with error bars
ax.errorbar(
    df["MWfrequencies"],
    df["SDmean"],
    yerr=df["err_SDmean"],
    fmt="o",
    color="black",
    markersize=5,
    capsize=2,
    elinewidth=1,
    markeredgewidth=0.5,
    zorder=2,
)

# Text box
textstr = r"Atom pumped in $F=1, m_F=0$"
props = dict(
    boxstyle="square,pad=0.4", facecolor="white", edgecolor="black", linewidth=1
)
ax.text(
    0.98,
    0.95,
    textstr,
    transform=ax.transAxes,
    fontsize=14,
    verticalalignment="top",
    horizontalalignment="right",
    bbox=props,
)

# Axes
ax.set_xlabel(r"MW frequency (kHz $-$ 6.834678 GHz)", fontsize=14)
ax.set_ylabel(r"Population in $F=2$", fontsize=14)
ax.tick_params(labelsize=13, direction="in", top=False, right=False)
ax.set_ylim(bottom=0, top=1.05)
ax.set_xlim(-200, 200)
ax.grid(False)

plt.tight_layout()
fig.savefig("MWSpec_F1mf0_thesis_style.pdf", dpi=300)
fig.savefig("MWSpec_F1mf0_thesis_style.png", dpi=300)
fig.savefig("MWSpec_F1mf0_thesis_style.svg")
# plt.show()
