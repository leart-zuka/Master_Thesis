import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update(
    {"font.size": 12, "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--"}
)


def plot_single(df, label, title, show_fig: bool = False):
    plt.figure(figsize=(8, 4.5))
    plt.plot(
        df["Time"],
        df["Ampl"],
        label=label,
        color="tab:blue",
        alpha=0.7,
        linewidth=2,
    )

    plt.xlabel("Time [ms]")
    plt.ylabel("Amplitude [V]")
    plt.title(title)
    plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(title, bbox_inches="tight", pad_inches=0.15)
    if show_fig:
        plt.show()


def plot_transmission(df_h, df_v, title):
    plt.figure(figsize=(8, 4.5))

    plt.plot(
        df_h["Time"],
        df_h["Ampl"],
        label="H pol.",
        color="tab:blue",
        alpha=0.7,
        linewidth=2,
    )

    plt.plot(
        df_v["Time"],
        df_v["Ampl"],
        label="V pol.",
        color="tab:orange",
        alpha=0.7,
        linewidth=2,
    )

    plt.xlabel("Time [ms]")
    plt.ylabel("Amplitude [V]")
    plt.title(title)
    plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(
        "coupling_for_same_input_power.svg", bbox_inches="tight", pad_inches=0.15
    )
    # plt.show()


df_h = pd.read_csv(
    "./C1H_transmission00001.csv",
    sep=",",
    skiprows=4,
)
plot_single(df_h, "H", "Tranmission Short Cavity for H", True)

df_v = pd.read_csv(
    "./C1V_transmission00001.csv",
    sep=",",
    skiprows=4,
)
plot_single(df_v, "V", "Tranmission Short Cavity for V", True)
df_a = pd.read_csv(
    "./C1A_transmission00001.csv",
    sep=",",
    skiprows=4,
)
plot_single(df_a, "A", "Tranmission Short Cavity for A", True)
df_d = pd.read_csv(
    "./C1D_transmission00001.csv",
    sep=",",
    skiprows=4,
)
plot_single(df_d, "D", "Tranmission Short Cavity for D", True)
df_r = pd.read_csv(
    "./C1R_transmission00001.csv",
    sep=",",
    skiprows=4,
)
plot_single(df_r, "R", "Tranmission Short Cavity for R", True)
df_l = pd.read_csv(
    "./C1L_transmission00001.csv",
    sep=",",
    skiprows=4,
)
plot_single(df_l, "L", "Tranmission Short Cavity for L", True)
