import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update(
    {"font.size": 12, "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--"}
)


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
    "./M1Polarization dependent loss for coupling between H (pink) and V (yellow)00000.csv",
    sep=",",
    skiprows=4,
)
df_v = pd.read_csv(
    "./M2Polarization dependent loss for coupling between H (pink) and V (yellow)00000.csv",
    sep=",",
    skiprows=4,
)

plot_transmission(
    df_h,
    df_v,
    "Comparisson of Transmission Signals for maximised V coupling for similar input powers",
)
