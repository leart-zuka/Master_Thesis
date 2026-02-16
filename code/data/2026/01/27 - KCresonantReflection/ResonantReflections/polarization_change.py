from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df_h = pd.read_csv("./H_start.csv", sep=";", skiprows=8)
df_v = pd.read_csv("./V_start.csv", sep=";", skiprows=8)
df_a = pd.read_csv("./A_start.csv", sep=";", skiprows=8)
df_d = pd.read_csv("./D_start.csv", sep=";", skiprows=8)
df_r = pd.read_csv("./R_start.csv", sep=";", skiprows=8)
df_l = pd.read_csv("./L_start.csv", sep=";", skiprows=8)


def plot_average_poincare_point(
    df_h,
    df_v,
    df_a,
    df_d,
    df_r,
    df_l,
    folder,
    mode: Literal["_reflection", "_transmission", ""],
    comparisson,
    s1_col=" Normalized s 1 ",
    s2_col=" Normalized s 2 ",
    s3_col=" Normalized s 3 ",
    show_plot=False,
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

    S1_a = df_a[s1_col].mean()
    S2_a = df_a[s2_col].mean()
    S3_a = df_a[s3_col].mean()

    S1_d = df_d[s1_col].mean()
    S2_d = df_d[s2_col].mean()
    S3_d = df_d[s3_col].mean()

    S1_r = df_r[s1_col].mean()
    S2_r = df_r[s2_col].mean()
    S3_r = df_r[s3_col].mean()

    S1_l = df_l[s1_col].mean()
    S2_l = df_l[s2_col].mean()
    S3_l = df_l[s3_col].mean()

    df_h_mode = pd.read_csv(f"./{folder}/H{mode}.csv", sep=";", skiprows=8)
    df_v_mode = pd.read_csv(f"./{folder}/V{mode}.csv", sep=";", skiprows=8)
    df_a_mode = pd.read_csv(f"./{folder}/A{mode}.csv", sep=";", skiprows=8)
    df_d_mode = pd.read_csv(f"./{folder}/D{mode}.csv", sep=";", skiprows=8)
    df_r_mode = pd.read_csv(f"./{folder}/R{mode}.csv", sep=";", skiprows=8)
    df_l_mode = pd.read_csv(f"./{folder}/L{mode}.csv", sep=";", skiprows=8)

    S1_h_compare = df_h_mode[s1_col].mean()
    S2_h_compare = df_h_mode[s2_col].mean()
    S3_h_compare = df_h_mode[s3_col].mean()

    S1_v_compare = df_v_mode[s1_col].mean()
    S2_v_compare = df_v_mode[s2_col].mean()
    S3_v_compare = df_v_mode[s3_col].mean()

    S1_a_compare = df_a_mode[s1_col].mean()
    S2_a_compare = df_a_mode[s2_col].mean()
    S3_a_compare = df_a_mode[s3_col].mean()

    S1_d_compare = df_d_mode[s1_col].mean()
    S2_d_compare = df_d_mode[s2_col].mean()
    S3_d_compare = df_d_mode[s3_col].mean()

    S1_r_compare = df_r_mode[s1_col].mean()
    S2_r_compare = df_r_mode[s2_col].mean()
    S3_r_compare = df_r_mode[s3_col].mean()

    S1_l_compare = df_l_mode[s1_col].mean()
    S2_l_compare = df_l_mode[s2_col].mean()
    S3_l_compare = df_l_mode[s3_col].mean()

    # ---- Sphere ----
    u = np.linspace(0, 2 * np.pi, 200)
    v = np.linspace(0, np.pi, 200)

    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    fig = plt.figure(figsize=(8, 6))
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
    #
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

    ax.scatter(S1_h, S2_h, S3_h, color="lightskyblue", s=80, label="Start H")
    ax.scatter(S1_v, S2_v, S3_v, color="lightcoral", s=80, label="Start V")
    ax.scatter(S1_d, S2_d, S3_d, color="plum", s=80, label="Start D")
    ax.scatter(S1_a, S2_a, S3_a, color="thistle", s=80, label="Start A")
    ax.scatter(S1_r, S2_r, S3_r, color="lightgreen", s=80, label="Start R")
    ax.scatter(S1_l, S2_l, S3_l, color="khaki", s=80, label="Start L")

    # ---- Comparison (AFTER) ----
    ax.scatter(
        S1_h_compare,
        S2_h_compare,
        S3_h_compare,
        color="dodgerblue",
        s=80,
        label=f"H after {mode}",
    )
    ax.scatter(
        S1_v_compare,
        S2_v_compare,
        S3_v_compare,
        color="firebrick",
        s=80,
        label=f"V after {mode}",
    )
    ax.scatter(
        S1_d_compare,
        S2_d_compare,
        S3_d_compare,
        color="purple",
        s=80,
        label=f"D after {mode}",
    )
    ax.scatter(
        S1_a_compare,
        S2_a_compare,
        S3_a_compare,
        color="indigo",
        s=80,
        label=f"A after {mode}",
    )
    ax.scatter(
        S1_r_compare,
        S2_r_compare,
        S3_r_compare,
        color="forestgreen",
        s=80,
        label=f"R after {mode}",
    )
    ax.scatter(
        S1_l_compare,
        S2_l_compare,
        S3_l_compare,
        color="goldenrod",
        s=80,
        label=f"L after {mode}",
    )
    # ---- Formatting ----
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()

    # plt.tight_layout()
    plt.title(comparisson)
    legend = ax.legend(loc="upper right", bbox_to_anchor=(1.05, 1), frameon=False)
    plt.savefig(f"{comparisson}.svg", bbox_inches="tight", pad_inches=0.2)
    if show_plot:
        plt.show()

    return S1_h, S2_h, S3_h, S1_v, S2_v, S3_v


plot_average_poincare_point(
    df_h,
    df_v,
    df_a,
    df_d,
    df_r,
    df_l,
    "polarimeterLTM",
    "",
    "Reflection_Resonant_Light",
    show_plot=True,
)
