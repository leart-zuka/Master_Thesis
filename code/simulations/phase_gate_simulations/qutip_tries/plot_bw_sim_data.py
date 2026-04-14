import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# --- Data ---
df = pd.read_csv("sim_data_bw_scan_scan.csv")
transmissions = df["Transmissions"].values
overlaps = df["Overlaps"].values

peak_idx = np.argmax(overlaps)
T_peak = transmissions[peak_idx]
overlap_peak = overlaps[peak_idx]

# --- Style ---
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
    }
)

fig, ax = plt.subplots(figsize=(5, 2.8))  # single-column thesis width

color = "#2872b0"

ax.plot(transmissions, overlaps, color=color, linewidth=1.4, zorder=2)
ax.scatter(transmissions, overlaps, color=color, s=12, zorder=3, linewidths=0)

# Peak annotation
ax.axvline(T_peak, color=color, linewidth=0.7, linestyle="--", alpha=0.5)
# ax.annotate(
#     rf"$T^* = {T_peak:.2f}$",
#     xy=(T_peak, overlap_peak),
#     xytext=(T_peak + 0.12, overlap_peak - 0.025),
#     fontsize=8,
#     color=color,
#     arrowprops=dict(arrowstyle="-", color=color, lw=0.7),
# )

ax.set_xlabel(r"V-mode transmission $T_\mathrm{V}$")
ax.set_ylabel("CNOT overlap")
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(0.42, 0.90)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.025))

fig.tight_layout()
fig.savefig("./plots/v_transmission_cnot_overlap.pdf")
fig.savefig("./plots/v_transmission_cnot_overlap.svg")
# plt.show()
