import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Define the sequence parameters (relative widths)
# MW pulse is longest, Detection is medium, Photon is short

labels = [
    "MW Pi Pulse\n" + r"$|0\rangle \leftrightarrow |1\rangle$",
    "Send\nPhoton",
    "State\nDetection",
]
# Pulse parameters (in microseconds)
durations = {"MW": 30, "Photon": 1, "SD": 7}

# Start positions
start_mw = 10
start_photon = start_mw + durations["MW"] + 10  # Gap of 10
start_sd = start_photon + durations["Photon"] + 10  # Gap of 10

fig, ax = plt.subplots(figsize=(10, 3))

# Define blocks: (label, start, duration)
blocks = [
    (
        "MW Pi Pulse\n" + r"$|0\rangle \leftrightarrow |1\rangle$",
        start_mw,
        durations["MW"],
        "#d3d3d3",
    ),
    ("Send Photon", start_photon, durations["Photon"], "#a8d5e2"),
    ("State \nDetection", start_sd, durations["SD"], "#f9dcc4"),
]

for label, start, duration, color in blocks:
    # Drawing the pulse block
    rect = patches.Rectangle(
        (start, 0.2),
        duration,
        0.6,
        linewidth=1.5,
        edgecolor="black",
        facecolor=color,
    )
    ax.add_patch(rect)

    # Label placement: inside for wide pulses, above for very narrow ones (like the 1us photon)
    if duration < 5:
        ax.text(
            start + duration / 2,
            0.85,
            label,
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    else:
        ax.text(
            start + duration / 2,
            0.5,
            label,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

# Axis configuration
ax.set_xlim(0, start_sd + durations["SD"] + 10)
ax.set_ylim(0, 1.2)
ax.set_xlabel("Time ($\mu$s)", fontsize=12, fontweight="bold")

# Hide the Y-axis and all spines except the bottom one
ax.get_yaxis().set_visible(False)
for spine in ["left", "top", "right"]:
    ax.spines[spine].set_visible(False)

# Add baseline
ax.axhline(0.2, color="black", linewidth=1.5)

plt.tight_layout()
plt.savefig("./plots/experimental_sequence.svg")
