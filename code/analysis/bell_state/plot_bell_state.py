import numpy as np
import matplotlib.pyplot as plt
from typing import Literal


PhotonBasis = Literal["Vpi", "RL"]


def _draw_wireframe_bar(ax, x: float, y: float, h: float, dx: float, dy: float) -> None:
    """Draw the 12 edges of a rectangular bar as blue lines."""
    # Corners of the box footprint (counter-clockwise)
    corners_bottom = [
        (x, y, 0),
        (x + dx, y, 0),
        (x + dx, y + dy, 0),
        (x, y + dy, 0),
    ]
    corners_top = [(cx, cy, h) for (cx, cy, _) in corners_bottom]

    edges = []
    # Bottom
    for i in range(4):
        edges.append((corners_bottom[i], corners_bottom[(i + 1) % 4]))
    # Top
    for i in range(4):
        edges.append((corners_top[i], corners_top[(i + 1) % 4]))
    # Verticals
    for i in range(4):
        edges.append((corners_bottom[i], corners_top[i]))

    for p0, p1 in edges:
        xs = [p0[0], p1[0]]
        ys = [p0[1], p1[1]]
        zs = [p0[2], p1[2]]
        ax.plot(xs, ys, zs, color="blue", linewidth=1.0, alpha=0.6)


def city_plot(
    rho: np.ndarray,
    *,
    photon_basis: PhotonBasis = "RL",
    target: np.ndarray | None = None,
    title: str = "",
    value: Literal["abs", "real", "imag"] = "abs",
    figsize: tuple[float, float] = (8, 6),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Render a 3D bar-chart city plot of a 4x4 density matrix.

    Parameters
    ----------
    rho : 4x4 complex density matrix in the {|0V>, |0pi>, |1V>, |1pi>} basis
          (as produced by density_matrix_from_stokes + project_to_physical).
    photon_basis : "Vpi" to plot as-is, "RL" to first rotate the photon to
          the R/L basis. Atom always stays in {|0>, |1>}.
    target : optional 4x4 ideal density matrix in the same basis; if given,
          drawn as open wireframe bars behind the solid data bars, Welte-style.
    value : which component of each matrix element to plot as bar height.
    """

    def _extract(m: np.ndarray) -> np.ndarray:
        if value == "abs":
            return np.abs(m)
        if value == "real":
            return m.real
        if value == "imag":
            return m.imag
        raise ValueError(f"Unknown value mode {value!r}")

    heights = _extract(rho)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    # Bar positions
    _xpos, _ypos = np.meshgrid(np.arange(4), np.arange(4), indexing="ij")
    xpos = _xpos.flatten()
    ypos = _ypos.flatten()
    zpos = np.zeros_like(xpos, dtype=float)
    dx = dy = 0.7 * np.ones_like(xpos, dtype=float)
    dz = heights.flatten()

    # Solid data bars
    ax.bar3d(
        xpos, ypos, zpos, dx, dy, dz, color="crimson", alpha=0.8, edgecolor="black"
    )

    # Optional wireframe target bars
    if target is not None:
        tgt_heights = _extract(target).flatten()
        for xi, yi, hi in zip(xpos, ypos, tgt_heights):
            if abs(hi) < 1e-9:
                continue
            # Draw a rectangular wireframe at the target height
            _draw_wireframe_bar(ax, xi, yi, hi, dx[0], dy[0])

    labels = [r"$|0R\rangle$", r"$|0L\rangle$", r"$|1R\rangle$", r"$|1L\rangle$"]
    ax.set_xticks(np.arange(4) + 0.35)
    ax.set_yticks(np.arange(4) + 0.35)
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.set_zlim(0, max(0.55, dz.max() * 1.1))
    ax.set_zlabel(
        {
            "abs": r"$|\rho_{ij}|$",
            "real": r"$\mathrm{Re}\,\rho_{ij}$",
            "imag": r"$\mathrm{Im}\,\rho_{ij}$",
        }[value]
    )
    if title:
        ax.set_title(title)

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig


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

city_plot(
    rho=rho_L,
    target=target_rho_L,
    title="test",
    value="abs",
    figsize=(10, 10),
)
