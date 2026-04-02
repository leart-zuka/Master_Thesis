import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rich import print
from rich.console import Console
from rich.table import Table

files = {
    "BW1": ("bw1_h.csv", "bw1_v.csv"),
    "BW2": ("bw2_h.csv", "bw2_v.csv"),
    "BW3": ("bw3_h.csv", "bw3_v.csv"),
    "BW4": ("bw4_h.csv", "bw4_v.csv"),
    "BW5": ("bw5_h.csv", "bw5_v.csv"),
    "BW6": ("bw6_h.csv", "bw6_v.csv"),
    "BW7": ("bw7_h.csv", "bw7_v.csv"),
    "BW8": ("bw8_h.csv", "bw8_v.csv"),
}

to_be_skipped_rows = 22


def convert_string_data_to_float(data: np.ndarray) -> np.ndarray:
    converted_data = np.zeros_like(data)
    for i, point in enumerate(data):
        converted_data[i] = float(point.replace(",", "."))
    return converted_data


baseline_df = pd.read_csv(
    "./data/baseline.csv", skiprows=to_be_skipped_rows, delimiter=";"
)
power_baseline = convert_string_data_to_float(baseline_df["Power (W)"].to_numpy())
avg_baseline = np.mean(power_baseline)
err_baseline = np.std(power_baseline)

print(f"Baseline Power is = {avg_baseline * 1e6:.2f} ± {err_baseline * 1e6:.2f} µW")

fig, axs = plt.subplots(4, 2, figsize=(12, 10))
axs = axs.flatten()

ax_index = 0

color_base = "C2"
color_H = "C3"
color_V = "C4"

style_baseline = "-"
style_bw = "--"
style_mean = ":"

transmission_h = []
transmission_h_err = []
transmission_v = []
transmission_v_err = []

for bw_name, (file_h, file_v) in files.items():
    df_h = pd.read_csv(f"./data/{file_h}", skiprows=to_be_skipped_rows, delimiter=";")
    ph = convert_string_data_to_float(df_h["Power (W)"].to_numpy())
    avg_h = np.mean(ph)
    err_h = np.std(ph)

    df_v = pd.read_csv(f"./data/{file_v}", skiprows=to_be_skipped_rows, delimiter=";")
    pv = convert_string_data_to_float(df_v["Power (W)"].to_numpy())
    avg_v = np.mean(pv)
    err_v = np.std(pv)

    # --- Transmittance (H) ---
    T_h = avg_h / avg_baseline
    dT_h = T_h * np.sqrt((err_h / avg_h) ** 2 + (err_baseline / avg_baseline) ** 2)

    transmission_h.append(T_h)
    transmission_h_err.append(dT_h)

    # --- Transmittance (V) ---
    T_v = avg_v / avg_baseline
    dT_v = T_v * np.sqrt((err_v / avg_v) ** 2 + (err_baseline / avg_baseline) ** 2)

    transmission_v.append(T_v)
    transmission_v_err.append(dT_v)

    console = Console()

    table = Table(show_header=False, box=None)

    console.print(f"[bold cyan]--- {bw_name} ---[/bold cyan]")

    t = Table(
        title=None,
        show_edge=False,
        show_header=False,
        pad_edge=False,
        box=None,
    )

    # H row
    t.add_row(
        "[bold blue]H[/bold blue]",
        f"{T_h * 100:.2f}% ± {dT_h * 100:.2f}%",
        f"{avg_h * 1000:.2f} ± {err_h * 1000:.2f} µW",
    )

    # V row
    t.add_row(
        "[bold red]V[/bold red]",
        f"{T_v * 100:.2f}% ± {dT_v * 100:.2f}%",
        f"{avg_v * 1000:.2f} ± {err_v * 1000:.2f} µW",
    )

    console.print(t)
    console.print("[dim]" + "─" * 42 + "[/dim]")

    # Plotting
    x_h = np.arange(len(ph))
    x_v = np.arange(len(pv))
    x_baseline = np.arange(len(power_baseline))

    ax = axs[ax_index]

    ax.plot(
        x_baseline,
        power_baseline,
        label="Baseline",
        color=color_base,
        linestyle=style_baseline,
        linewidth=2,
        alpha=0.7,
    )

    ax.plot(
        x_h,
        ph,
        label=f"{bw_name} – H",
        color=color_H,
        linestyle=style_bw,
        linewidth=2,
    )

    ax.plot(
        x_v,
        pv,
        label=f"{bw_name} – V",
        color=color_V,
        linestyle=style_bw,
        linewidth=2,
    )

    ax.axhline(avg_baseline, color=color_base, linestyle=style_mean, alpha=0.5)
    ax.axhline(avg_h, color=color_H, linestyle=style_mean, alpha=0.9)
    ax.axhline(avg_v, color=color_V, linestyle=style_mean, alpha=0.9)

    ax.set_title(f"{bw_name}: Power vs Time")
    ax.set_xlabel("Sample index")
    ax.set_ylabel("Power [mW]")
    ax.legend(loc="center right")
    ax.grid(True)

    ax_index += 1

print(f"{np.mean(transmission_h)} +/- {np.mean(transmission_h_err)}")
print(f"{np.mean(transmission_v)} +/- {np.mean(transmission_v_err)}")

plt.tight_layout()
plt.title("Characterization of single Brewster Windows")
plt.savefig("single_brewster_window_characterization.pdf")
plt.savefig("single_brewster_window_characterization.svg")
plt.show()
