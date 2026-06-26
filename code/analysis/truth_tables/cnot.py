import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mticker

# ── RAW DATA ──
cnot_truth_table = np.array(
    [
        [0.85459184, 0.08291457, 0.0, 0.00295858],
        [0.08418367, 0.8241206, 0.0, 0.00887574],
        [0.04336735, 0.04271357, 0.06086957, 0.90828402],
        [0.01785714, 0.05025126, 0.93913043, 0.07988166],
    ]
)
cnot_truth_table_err = np.array(
    [
        [0.01780454, 0.01382225, 0.0, 0.0029542],
        [0.0140241, 0.01908365, 0.0, 0.00510162],
        [0.01028753, 0.0101359, 0.01576519, 0.01569911],
        [0.00668883, 0.01095056, 0.01576519, 0.01474643],
    ]
)

cnot_labels = [
    r"$|1\rangle \otimes R$",
    r"$|1\rangle \otimes L$",
    r"$|0\rangle \otimes R$",
    r"$|0\rangle \otimes L$",
]

cnot_ideal = np.array(
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float
)

bar_color = "#c97a2a"
bar_edge = "#1a1a1a"
ideal_face = "#8bbae0"
ideal_alpha = 0.20
ideal_edge = "#6a9ec4"

# ── FONT SIZES ──
TICK_FS = 12
LABEL_FS = 14
PANEL_FS = 16
TEXT_FS = 11  # in-cell numbers on heatmap

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.labelsize": LABEL_FS,
        "xtick.labelsize": TICK_FS,
        "ytick.labelsize": TICK_FS,
    }
)


def plot_3d_bars(ax, matrix, err, labels, ideal):
    n = 4
    xpos, ypos = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    xpos_f = xpos.flatten()
    ypos_f = ypos.flatten()
    dx = dy = 0.48
    dz = matrix.flatten()

    for i in range(n):
        for j in range(n):
            if ideal[i, j] > 0.5:
                ax.bar3d(
                    i + (1 - dx) / 2,
                    j + (1 - dy) / 2,
                    0,
                    dx,
                    dy,
                    1.0,
                    color=ideal_face,
                    alpha=ideal_alpha,
                    edgecolor=ideal_edge,
                    linewidth=0.6,
                    zsort="min",
                )

    ax.bar3d(
        xpos_f + (1 - dx) / 2,
        ypos_f + (1 - dy) / 2,
        np.zeros(n * n),
        dx,
        dy,
        dz,
        color=bar_color,
        edgecolor=bar_edge,
        alpha=0.88,
        linewidth=0.4,
        zsort="average",
    )

    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            e = err[i, j]
            if e > 0:
                cx, cy = i + 0.5, j + 0.5
                ax.plot(
                    [cx, cx],
                    [cy, cy],
                    [val - e, val + e],
                    color="k",
                    linewidth=0.7,
                    zorder=10,
                )
                cap = 0.08
                for z in [val + e, val - e]:
                    ax.plot(
                        [cx - cap, cx + cap],
                        [cy, cy],
                        [z, z],
                        color="k",
                        linewidth=0.5,
                        zorder=10,
                    )

    tick_pos = np.arange(n) + 0.5
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(
        labels, rotation=-20, ha="center", fontsize=TICK_FS, rotation_mode="anchor"
    )
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(
        labels, rotation=20, ha="center", fontsize=TICK_FS, rotation_mode="anchor"
    )
    ax.set_zlim(0, 1.0)
    ax.set_zticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.zaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax.tick_params(axis="z", labelsize=TICK_FS, pad=1)
    ax.tick_params(axis="x", pad=14)
    ax.tick_params(axis="y", pad=14)

    ax.set_xlabel("Output states", labelpad=32, fontsize=LABEL_FS)
    ax.set_ylabel("Input states", labelpad=32, fontsize=LABEL_FS)
    ax.set_zlabel("Normalized probability", labelpad=20, fontsize=LABEL_FS)

    ax.view_init(elev=25, azim=-50)
    ax.set_box_aspect(None, zoom=0.78)

    ax.xaxis.pane.set_facecolor((0.95, 0.95, 0.95, 0.5))
    ax.yaxis.pane.set_facecolor((0.90, 0.90, 0.90, 0.5))
    ax.zaxis.pane.set_facecolor((0.93, 0.93, 0.93, 0.5))
    ax.xaxis.pane.set_edgecolor((0.90, 0.90, 0.90, 1))
    ax.yaxis.pane.set_edgecolor((0.85, 0.85, 0.85, 1))
    ax.zaxis.pane.set_edgecolor((0.88, 0.88, 0.88, 1))
    ax.grid(True, alpha=0.3)


def plot_heatmap(ax, matrix, err, labels):
    n = 4
    im = ax.imshow(
        matrix,
        cmap="YlOrBr",
        vmin=0,
        vmax=1,
        aspect="equal",
        origin="upper",
        interpolation="nearest",
    )
    for i in range(n):
        for j in range(n):
            val = matrix[i, j]
            e = err[i, j]
            color = "white" if val > 0.45 else "#333333"
            fontw = "bold" if val > 0.5 else "normal"
            txt = f"{val:.2f}\n±{e:.2f}"
            ax.text(
                j,
                i,
                txt,
                ha="center",
                va="center",
                fontsize=TEXT_FS,
                fontweight=fontw,
                color=color,
                linespacing=1.3,
            )
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(labels, fontsize=TICK_FS)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(labels, fontsize=TICK_FS)
    ax.set_xlabel("Input states", fontsize=LABEL_FS)
    ax.set_ylabel("Output states", fontsize=LABEL_FS)
    for edge in np.arange(-0.5, n, 1):
        ax.axhline(edge, color="white", linewidth=2)
        ax.axvline(edge, color="white", linewidth=2)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#666666")
        spine.set_linewidth(0.8)
    return im


# ── FIGURE (a): 3D bars ──
fig_a = plt.figure(figsize=(10, 8))
ax_a = fig_a.add_subplot(111, projection="3d")
plot_3d_bars(ax_a, cnot_truth_table, cnot_truth_table_err, cnot_labels, cnot_ideal)
ax_a.text2D(
    0.05,
    0.95,
    "(a)",
    transform=ax_a.transAxes,
    fontsize=PANEL_FS,
    fontweight="bold",
    va="top",
)
fig_a.savefig("cnot_truth_table_3d.svg", bbox_inches="tight")
fig_a.savefig("cnot_truth_table_3d.png", bbox_inches="tight", dpi=250)
plt.close(fig_a)

# ── FIGURE (b): heatmap ──
fig_b, ax_b = plt.subplots(figsize=(7, 6))
im_b = plot_heatmap(ax_b, cnot_truth_table, cnot_truth_table_err, cnot_labels)
ax_b.text(
    -0.15,
    1.12,
    "(b)",
    transform=ax_b.transAxes,
    fontsize=PANEL_FS,
    fontweight="bold",
    va="top",
)
fig_b.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.15)
fig_b.savefig("cnot_truth_table_heatmap.svg", bbox_inches="tight")
fig_b.savefig("cnot_truth_table_heatmap.png", bbox_inches="tight", dpi=250)
plt.close(fig_b)

print("CNOT figures done")
