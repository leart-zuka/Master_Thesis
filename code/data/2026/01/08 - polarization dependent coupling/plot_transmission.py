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

    plt.savefig(f"{title}.svg", bbox_inches="tight", pad_inches=0.15)
    plt.show()


df_h = pd.read_csv(
    "./C1H_transmission_maximised_coupling_h00000.csv", sep=",", skiprows=4
)
df_v = pd.read_csv(
    "./C1V_transmission_maximised_coupling_h00000.csv", sep=",", skiprows=4
)

plot_transmission(
    df_h, df_v, "Transmission spectrum (coupling maximised for H polarization)"
)


df_h = pd.read_csv(
    "./C1H_transmission_maximised_coupling_v00000.csv", sep=",", skiprows=4
)
df_v = pd.read_csv(
    "./C1V_transmission_maximised_coupling_v00000.csv", sep=",", skiprows=4
)

plot_transmission(
    df_h, df_v, "Transmission spectrum (coupling maximised for V polarization)"
)
