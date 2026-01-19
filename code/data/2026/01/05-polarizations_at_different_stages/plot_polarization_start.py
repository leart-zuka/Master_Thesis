import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


df_h = pd.read_csv("./data/H_before_bw.csv", sep=";", skiprows=8)
df_v = pd.read_csv("./data/V_before_bw.csv", sep=";", skiprows=8)


def plot_average_poincare_point(
    df_h,
    df_v,
    df_h_compare,
    df_v_compare,
    comparisson,
    s1_col=" Normalized s 1 ",
    s2_col=" Normalized s 2 ",
    s3_col=" Normalized s 3 ",
):
    """
    Compute the average normalized Stokes vector and plot it on the Poincaré sphere.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing normalized Stokes parameters.
    s1_col, s2_col, s3_col : str
        Column names for normalized S1, S2, S3.
    """

    # ---- compute averages ----
    S1_h = df_h[s1_col].mean()
    S2_h = df_h[s2_col].mean()
    S3_h = df_h[s3_col].mean()

    S1_v = df_v[s1_col].mean()
    S2_v = df_v[s2_col].mean()
    S3_v = df_v[s3_col].mean()

    S1_h_compare = df_h_compare[s1_col].mean()
    S2_h_compare = df_h_compare[s2_col].mean()
    S3_h_compare = df_h_compare[s3_col].mean()

    S1_v_compare = df_v_compare[s1_col].mean()
    S2_v_compare = df_v_compare[s2_col].mean()
    S3_v_compare = df_v_compare[s3_col].mean()

    # ---- Sphere ----
    u = np.linspace(0, 2 * np.pi, 200)
    v = np.linspace(0, np.pi, 200)

    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Main sphere
    ax.plot_surface(x, y, z, color="#1f4bd8", alpha=0.08, linewidth=0)

    # ---- Reference great circles ----
    t = np.linspace(0, 2 * np.pi, 400)

    # S1–S2 plane (equator)
    ax.plot(np.cos(t), np.sin(t), 0 * t, color="gray", lw=1)

    # S1–S3 plane
    ax.plot(np.cos(t), 0 * t, np.sin(t), color="gray", lw=1)

    # S2–S3 plane
    ax.plot(0 * t, np.cos(t), np.sin(t), color="gray", lw=1)

    # ---- Axes ----
    ax.plot([-1, 1], [0, 0], [0, 0], color="black", lw=1)  # S1
    ax.plot([0, 0], [-1, 1], [0, 0], color="black", lw=1)  # S2
    ax.plot([0, 0], [0, 0], [-1, 1], color="black", lw=1)  # S3

    # ---- Canonical polarization points ----
    points = {
        "H": (1, 0, 0),
        "V": (-1, 0, 0),
        "D": (0, 1, 0),
        "A": (0, -1, 0),
        "R": (0, 0, 1),
        "L": (0, 0, -1),
    }

    colors = {
        "H": "green",
        "V": "green",
        "D": "orange",
        "A": "orange",
        "R": "red",
        "L": "red",
    }

    for label, (x0, y0, z0) in points.items():
        # ax.scatter(x0, y0, z0, color=colors[label], s=40)
        ax.text(
            1.1 * x0,
            1.1 * y0,
            1.1 * z0,
            label,
            color=colors[label],
            fontsize=12,
            weight="bold",
        )

    # ---- Average polarization point ----
    ax.scatter(S1_h, S2_h, S3_h, color="lightgreen", s=80, label="Start H")
    ax.scatter(S1_v, S2_v, S3_v, color="lightcoral", s=80, label="Start V")

    ax.scatter(
        S1_h_compare,
        S2_h_compare,
        S3_h_compare,
        color="darkgreen",
        s=80,
        label=f"H after {comparisson}",
    )
    ax.scatter(
        S1_v_compare,
        S2_v_compare,
        S3_v_compare,
        color="darkred",
        s=80,
        label=f"V after {comparisson}",
    )

    # ---- Formatting ----
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()

    plt.tight_layout()
    plt.title(comparisson)
    plt.legend()
    plt.show()

    return S1_h, S2_h, S3_h, S1_v, S2_v, S3_v


comparisson = "transmission through 1st dichroic"
df_h_compare = pd.read_csv(
    "./data/H_transmission_1st_dichroic.csv", sep=";", skiprows=8
)
df_v_compare = pd.read_csv(
    "./data/V_transmission_1st_dichroic.csv", sep=";", skiprows=8
)

plot_average_poincare_point(df_h, df_v, df_h_compare, df_v_compare, comparisson)
