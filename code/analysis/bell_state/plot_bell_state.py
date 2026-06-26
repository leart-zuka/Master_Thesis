import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
from typing import Literal


PhotonBasis = Literal["Vpi", "RL"]

# ------------------------------------------------------------------
# Font sizes
# ------------------------------------------------------------------
TICK_FS = 18  # axis tick labels (ket labels)
LABEL_FS = 22  # axis titles and z-label
PANEL_FS = 24  # (a) / (b) panel letters
CBAR_FS = 20  # colorbar label
TEXT_FS = 16  # in-cell numbers on heatmap


def _draw_wireframe_bar(ax, x: float, y: float, h: float, dx: float, dy: float) -> None:
    """Draw the 12 edges of a rectangular bar as blue lines."""
    corners_bottom = [
        (x, y, 0),
        (x + dx, y, 0),
        (x + dx, y + dy, 0),
        (x, y + dy, 0),
    ]
    corners_top = [(cx, cy, h) for (cx, cy, _) in corners_bottom]

    edges = []
    for i in range(4):
        edges.append((corners_bottom[i], corners_bottom[(i + 1) % 4]))
    for i in range(4):
        edges.append((corners_top[i], corners_top[(i + 1) % 4]))
    for i in range(4):
        edges.append((corners_bottom[i], corners_top[i]))

    for p0, p1 in edges:
        xs = [p0[0], p1[0]]
        ys = [p0[1], p1[1]]
        zs = [p0[2], p1[2]]
        ax.plot(xs, ys, zs, color="blue", linewidth=1.0, alpha=0.6)


def _extract(m: np.ndarray, value: str) -> np.ndarray:
    if value == "abs":
        return np.abs(m)
    if value == "abs2":
        return np.abs(m) ** 2
    if value == "real":
        return m.real
    if value == "imag":
        return m.imag
    raise ValueError(f"Unknown value mode {value!r}")


def city_plot_and_heatmap(
    rho: np.ndarray,
    *,
    target: np.ndarray | None = None,
    title: str = "",
    value: Literal["abs", "abs2", "real", "imag"] = "abs2",
    figsize: tuple[float, float] = (14, 7),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Render a side-by-side figure: 3D city plot on the left, heatmap on the right.
    Also saves each subplot as a separate file when save_path is given.
    """

    _value_label = {
        "abs": r"$|\rho_{ij}|$",
        "abs2": r"$|\rho_{ij}|^2$",
        "real": r"$\mathrm{Re}\,\rho_{ij}$",
        "imag": r"$\mathrm{Im}\,\rho_{ij}$",
    }[value]

    heights = _extract(rho, value)
    labels = [r"$|0R\rangle$", r"$|0L\rangle$", r"$|1R\rangle$", r"$|1L\rangle$"]

    # ---------------------------------------------------------------
    # Helper: build the 3D city axis in any figure
    # ---------------------------------------------------------------
    def _make_3d(ax3d):
        _xpos, _ypos = np.meshgrid(np.arange(4), np.arange(4), indexing="ij")
        xpos = _xpos.flatten()
        ypos = _ypos.flatten()
        zpos = np.zeros_like(xpos, dtype=float)
        dx = dy = 0.7 * np.ones_like(xpos, dtype=float)
        dz = heights.flatten()

        ax3d.bar3d(
            xpos, ypos, zpos, dx, dy, dz, color="crimson", alpha=0.8, edgecolor="black"
        )

        if target is not None:
            tgt_heights = _extract(target, value).flatten()
            for xi, yi, hi in zip(xpos, ypos, tgt_heights):
                if abs(hi) < 1e-9:
                    continue
                _draw_wireframe_bar(ax3d, xi, yi, hi, dx[0], dy[0])

        ax3d.set_xticks(np.arange(4) + 0.35)
        ax3d.set_yticks(np.arange(4) + 0.35)
        ax3d.set_xticklabels(labels, fontsize=TICK_FS)
        ax3d.set_yticklabels(labels, fontsize=TICK_FS)
        ax3d.set_zlim(0, 0.50)
        ax3d.set_zlabel(_value_label, labelpad=15, rotation=90, fontsize=LABEL_FS)
        ax3d.tick_params(axis="z", labelsize=TICK_FS)
        ax3d.view_init(elev=15, azim=-30)

    # Helper: build the 2D heatmap axis
    def _make_2d(ax2d):
        vmax = heights.max()
        if target is not None:
            vmax = max(vmax, _extract(target, value).max())
        vmax = max(vmax, 1e-9)

        im = ax2d.imshow(
            heights,
            cmap="Reds",
            norm=Normalize(vmin=0, vmax=vmax),
            aspect="equal",
            origin="upper",
        )

        threshold = vmax * 0.5
        for i in range(4):
            for j in range(4):
                val = heights[i, j]
                color = "white" if val > threshold else "black"
                ax2d.text(
                    j,
                    i,
                    f"{val:.3f}",
                    ha="center",
                    va="center",
                    color=color,
                    fontsize=TEXT_FS,
                )

        ax2d.set_xticks(np.arange(4))
        ax2d.set_yticks(np.arange(4))
        ax2d.set_xticklabels(labels, fontsize=TICK_FS)
        ax2d.set_yticklabels(labels, fontsize=TICK_FS)
        ax2d.set_xlabel("Column index", fontsize=LABEL_FS)
        ax2d.set_ylabel("Row index", fontsize=LABEL_FS)
        return im, vmax

    # ---------------------------------------------------------------
    # Combined figure
    # ---------------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 2, width_ratios=[1.4, 1.0], figure=fig)

    ax3d = fig.add_subplot(gs[0, 0], projection="3d")
    _make_3d(ax3d)

    ax3d_bbox = ax3d.get_position()
    fig.text(
        ax3d_bbox.x0 + 0.01,
        ax3d_bbox.y1 - 0.03,
        "(a)",
        fontsize=PANEL_FS,
        fontweight="bold",
        ha="left",
        va="top",
    )

    ax2d = fig.add_subplot(gs[0, 1])
    im, vmax = _make_2d(ax2d)

    ax2d.text(
        -0.15,
        1.05,
        "(b)",
        transform=ax2d.transAxes,
        fontsize=PANEL_FS,
        fontweight="bold",
        ha="left",
        va="top",
    )

    cbar = fig.colorbar(im, ax=ax2d, fraction=0.046, pad=0.04)
    cbar.set_label(_value_label, fontsize=CBAR_FS)
    cbar.ax.tick_params(labelsize=TICK_FS)

    if title:
        fig.suptitle(title, y=0.98, fontsize=13)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.08, wspace=0.30)

    if save_path:
        for ext in ("png", "pdf", "svg"):
            fig.savefig(f"./{save_path}.{ext}", dpi=200, bbox_inches="tight")

        # -----------------------------------------------------------
        # Subplot (a): 3D city plot — standalone
        # -----------------------------------------------------------
        fig_a = plt.figure(figsize=(8, 7))
        ax3d_a = fig_a.add_subplot(111, projection="3d")
        _make_3d(ax3d_a)
        fig_a.text(
            0.02,
            0.96,
            "(a)",
            fontsize=PANEL_FS,
            fontweight="bold",
            ha="left",
            va="top",
        )
        if title:
            fig_a.suptitle(title, y=0.98, fontsize=13)
        fig_a.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.08)
        for ext in ("png", "pdf", "svg"):
            fig_a.savefig(f"./{save_path}_3d.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig_a)

        # -----------------------------------------------------------
        # Subplot (b): heatmap — standalone
        # -----------------------------------------------------------
        fig_b, ax2d_b = plt.subplots(figsize=(7, 6))
        im_b, _ = _make_2d(ax2d_b)
        ax2d_b.text(
            -0.15,
            1.05,
            "(b)",
            transform=ax2d_b.transAxes,
            fontsize=PANEL_FS,
            fontweight="bold",
            ha="left",
            va="top",
        )
        cbar_b = fig_b.colorbar(im_b, ax=ax2d_b, fraction=0.046, pad=0.04)
        cbar_b.set_label(_value_label, fontsize=CBAR_FS)
        cbar_b.ax.tick_params(labelsize=TICK_FS)
        if title:
            fig_b.suptitle(title, y=0.98, fontsize=13)
        fig_b.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.12)
        for ext in ("png", "pdf", "svg"):
            fig_b.savefig(f"./{save_path}_heatmap.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig_b)

    return fig


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

target_rho_L = np.array(
    [
        [0.5 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, -0.5 + 0.0j],
        [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j],
        [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j],
        [-0.5 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 0.5 + 0.0j],
    ]
)

rho_L = np.array(
    [
        [
            0.35190793 + 0.0j,
            0.01840359 + 0.00322469j,
            0.00624527 - 0.02268682j,
            -0.28188186 - 0.04753902j,
        ],
        [
            0.01840359 - 0.00322469j,
            0.11207603 + 0.0j,
            0.02305752 + 0.11666036j,
            -0.03840008 + 0.08026703j,
        ],
        [
            0.00624527 + 0.02268682j,
            0.02305752 - 0.11666036j,
            0.13272999 + 0.0j,
            0.09010545 + 0.02156167j,
        ],
        [
            -0.28188186 + 0.04753902j,
            -0.03840008 - 0.08026703j,
            0.09010545 - 0.02156167j,
            0.40328605 + 0.0j,
        ],
    ]
)


if __name__ == "__main__":
    city_plot_and_heatmap(
        rho=rho_L,
        target=target_rho_L,
        title="",
        value="abs",
        figsize=(14, 7),
        save_path="bell_state_L",
    )
