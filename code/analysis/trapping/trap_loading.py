import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

df = pd.read_csv("atom_6_loading.csv")

time_ms = df["time_s"] * 1000
counts = df["counts"]

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 11,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.major.size": 5,
        "ytick.minor.size": 3,
        "xtick.top": True,
        "ytick.right": True,
    }
)

fig, ax = plt.subplots(figsize=(6, 3.2), dpi=200)

ax.bar(
    time_ms,
    counts,
    width=20,
    align="edge",
    color="#3a6ea5",
    edgecolor="#3a6ea5",
    linewidth=0.3,
)

# --- Phase dividers (adjust these times in ms) ---
dividers = [500, 1000, 1300]

for x in dividers:
    ax.axvline(x, color="grey", linewidth=0.8, linestyle="-", zorder=3)

# --- Circled labels between dividers ---
boundaries = [0] + dividers + [time_ms.max() + 20]
label_y = counts.max() * 0.92
for i in range(len(boundaries) - 1):
    mid = (boundaries[i] + boundaries[i + 1]) / 2
    ax.text(
        mid,
        0.95,
        str(i + 1),
        ha="center",
        va="center",
        fontsize=7,
        color="grey",
        transform=ax.get_xaxis_transform(),
        bbox=dict(
            boxstyle="circle,pad=0.3",
            facecolor="white",
            edgecolor="grey",
            linewidth=0.8,
        ),
    )

ax.set_xlabel("Time (ms)")
ax.set_ylabel("Fluorescence counts / 20 ms")

ax.xaxis.set_major_locator(ticker.MultipleLocator(500))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1000))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(500))

ax.set_xlim(0, time_ms.max() + 20)
ax.set_ylim(0, None)

fig.tight_layout()
fig.savefig("atom_6_loading.svg", bbox_inches="tight")
fig.savefig("atom_6_loading.pdf", bbox_inches="tight")
# plt.show()
