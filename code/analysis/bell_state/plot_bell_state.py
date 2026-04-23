import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
from typing import Literal


PhotonBasis = Literal["Vpi", "RL"]


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

    Parameters
    ----------
    rho : 4x4 complex density matrix in the {|0R>, |0L>, |1R>, |1L>} basis
          (assumes it's already been rotated into the R/L display basis).
    target : optional 4x4 ideal density matrix in the same basis; drawn as
          open wireframe bars behind the solid data bars on the city plot.
    value : which component of each matrix element to plot.
            "abs"  -> |rho_ij|
            "abs2" -> |rho_ij|^2
            "real" -> Re(rho_ij)
            "imag" -> Im(rho_ij)
    """

    def _extract(m: np.ndarray) -> np.ndarray:
        if value == "abs":
            return np.abs(m)
        if value == "abs2":
            return np.abs(m) ** 2
        if value == "real":
            return m.real
        if value == "imag":
            return m.imag
        raise ValueError(f"Unknown value mode {value!r}")

    _value_label = {
        "abs": r"$|\rho_{ij}|$",
        "abs2": r"$|\rho_{ij}|^2$",
        "real": r"$\mathrm{Re}\,\rho_{ij}$",
        "imag": r"$\mathrm{Im}\,\rho_{ij}$",
    }[value]

    heights = _extract(rho)
    labels = [r"$|0R\rangle$", r"$|0L\rangle$", r"$|1R\rangle$", r"$|1L\rangle$"]

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 2, width_ratios=[1.4, 1.0], figure=fig)

    # ------------------------------------------------------------------
    # Left: 3D city plot
    # ------------------------------------------------------------------
    ax3d = fig.add_subplot(gs[0, 0], projection="3d")

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
        tgt_heights = _extract(target).flatten()
        for xi, yi, hi in zip(xpos, ypos, tgt_heights):
            if abs(hi) < 1e-9:
                continue
            _draw_wireframe_bar(ax3d, xi, yi, hi, dx[0], dy[0])

    ax3d.set_xticks(np.arange(4) + 0.35)
    ax3d.set_yticks(np.arange(4) + 0.35)
    ax3d.set_xticklabels(labels)
    ax3d.set_yticklabels(labels)
    # z-axis upper limit should be set by the larger of data or target
    if target is not None:
        z_max = max(dz.max(), _extract(target).max()) * 1.1
    else:
        z_max = dz.max() * 1.1
    ax3d.set_zlim(0, max(0.50, z_max))
    ax3d.set_zlabel(_value_label, labelpad=15, rotation=90)
    ax3d.view_init(elev=15, azim=-30)

    # ------------------------------------------------------------------
    # Right: 2D heatmap
    # ------------------------------------------------------------------
    ax2d = fig.add_subplot(gs[0, 1])

    # Use same color scale as city plot for visual consistency
    vmax = heights.max()
    if target is not None:
        vmax = max(vmax, _extract(target).max())
    vmax = max(vmax, 1e-9)  # guard against all-zero case

    im = ax2d.imshow(
        heights,
        cmap="Reds",
        norm=Normalize(vmin=0, vmax=vmax),
        aspect="equal",
        origin="upper",
    )

    # Overlay numeric values in each cell for quantitative readability
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
                fontsize=10,
            )

    ax2d.set_xticks(np.arange(4))
    ax2d.set_yticks(np.arange(4))
    ax2d.set_xticklabels(labels)
    ax2d.set_yticklabels(labels)
    ax2d.set_xlabel("Column index")
    ax2d.set_ylabel("Row index")

    # Colorbar
    cbar = fig.colorbar(im, ax=ax2d, fraction=0.046, pad=0.04)
    cbar.set_label(_value_label)

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------
    if title:
        fig.suptitle(title, y=0.98, fontsize=13)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.08, wspace=0.30)

    if save_path:
        fig.savefig(f"./{save_path}.png", dpi=200, bbox_inches="tight")
        fig.savefig(f"./{save_path}.pdf", bbox_inches="tight")
        fig.savefig(f"./{save_path}.svg", bbox_inches="tight")

    return fig


# ---------------------------------------------------------------------------
# Example data (same as your script)
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
    # plt.show()
