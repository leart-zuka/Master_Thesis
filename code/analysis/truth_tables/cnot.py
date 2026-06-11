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

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.labelsize": 9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
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
    ax.set_xticklabels(labels, rotation=-25, ha="left", fontsize=7)
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(labels, rotation=25, ha="right", fontsize=7)
    ax.set_zlim(0, 1.0)
    ax.set_zticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.zaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax.tick_params(axis="z", labelsize=7, pad=1)

    ax.set_xlabel("Output states", labelpad=10, fontsize=8)
    ax.set_ylabel("Input states", labelpad=10, fontsize=8)
    ax.set_zlabel("Normalized probability", labelpad=2, fontsize=8)

    ax.view_init(elev=25, azim=-50)
    ax.set_box_aspect(None, zoom=0.85)

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
                fontsize=7.5,
                fontweight=fontw,
                color=color,
                linespacing=1.3,
            )
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Input states", fontsize=9)
    ax.set_ylabel("Output states", fontsize=9)
    for edge in np.arange(-0.5, n, 1):
        ax.axhline(edge, color="white", linewidth=2)
        ax.axvline(edge, color="white", linewidth=2)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#666666")
        spine.set_linewidth(0.8)
    return im


# ── BUILD FIGURE ──
fig = plt.figure(figsize=(11, 5))
gs = fig.add_gridspec(
    1,
    2,
    hspace=0.3,
    wspace=0.30,
    bottom=0.10,
    top=0.95,
    left=0.05,
    right=0.95,
)

ax1 = fig.add_subplot(gs[0, 0], projection="3d")
plot_3d_bars(ax1, cnot_truth_table, cnot_truth_table_err, cnot_labels, cnot_ideal)

ax2 = fig.add_subplot(gs[0, 1])
im = plot_heatmap(ax2, cnot_truth_table, cnot_truth_table_err, cnot_labels)

for ax_obj, label in zip([ax1, ax2], ["(a)", "(b)"]):
    if hasattr(ax_obj, "zaxis"):
        ax_obj.text2D(
            0.05,
            0.95,
            label,
            transform=ax_obj.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )
    else:
        ax_obj.text(
            -0.15,
            1.12,
            label,
            transform=ax_obj.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )

# cbar_ax = fig.add_axes([0.15, 0.015, 0.7, 0.012])
# cbar = fig.colorbar(im, cax=cbar_ax, pad=0.3, orientation="horizontal")
# cbar.set_label("Normalized probability", fontsize=9)
# cbar.ax.tick_params(labelsize=8)

plt.savefig("cnot_truth_table.svg", bbox_inches="tight")
plt.savefig("cnot_truth_table.png", bbox_inches="tight", dpi=250)
print("CNOT figure done")
