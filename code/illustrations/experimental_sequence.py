import matplotlib.pyplot as plt
import matplotlib.patches as patches

# 1. Configuration
real_durations = {"OP": 600, "MW": 30, "Photon": 1, "SD": 7}
# Visual widths adjusted to fit text comfortably
visual_durations = {"OP": 45, "MW": 35, "Photon": 15, "SD": 22}
gap = 4

# Calculate positions
start_op = 2
start_mw = start_op + visual_durations["OP"] + gap
start_photon = start_mw + visual_durations["MW"] + gap
start_sd = start_photon + visual_durations["Photon"] + gap

fig, ax = plt.subplots(figsize=(10, 2.5))

# 2. Define blocks: (Display Name, Timing, Start, Visual Width, Color)
blocks = [
    (
        "OP + Cooling",
        f"{real_durations['OP']} $\mu$s",
        start_op,
        visual_durations["OP"],
        "#c7e9c0",
    ),
    (
        "MW $\pi$ Pulse" + "\n" + r"$|0\rangle \leftrightarrow |1\rangle$",
        f"{real_durations['MW']} $\mu$s",
        start_mw,
        visual_durations["MW"],
        "#d3d3d3",
    ),
    (
        "Photon",
        f"{real_durations['Photon']} $\mu$s",
        start_photon,
        visual_durations["Photon"],
        "#a8d5e2",
    ),
    (
        "State Detection",
        f"{real_durations['SD']} $\mu$s",
        start_sd,
        visual_durations["SD"],
        "#f9dcc4",
    ),
]

for label, time_val, start, v_dur, color in blocks:
    # Draw the pulse block
    rect = patches.Rectangle(
        (start, 0.3),
        v_dur,
        0.5,
        linewidth=1.5,
        edgecolor="black",
        facecolor=color,
        zorder=2,
    )
    ax.add_patch(rect)

    # Place Name and Timing inside the block
    ax.text(
        start + v_dur / 2,
        0.55,
        f"{label}\n({time_val})",
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
    )

# 3. Clean Axis Styling
ax.set_xlim(0, start_sd + visual_durations["SD"] + 5)
ax.set_ylim(0, 1.1)

# Remove all ticks and numerical values
ax.set_xticks([])
ax.set_yticks([])

# Set "Time" as the simple label
ax.set_xlabel("Time", fontsize=12, fontweight="bold")

# Hide all border spines except the bottom one
for spine in ["left", "top", "right"]:
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_linewidth(1.5)

# Add baseline
ax.axhline(0.3, color="black", linewidth=1.5, zorder=1)

plt.tight_layout()
plt.savefig("./plots/experimental_sequence.svg")
