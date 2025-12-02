import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rich.console import Console
from rich.table import Table

baseline_h_file = "baseline_h.csv"
baseline_v_file = "baseline_v.csv"

plt.style.use("seaborn-v0_8")

color_H_base = "C2"
color_V_base = "C1"
color_H = "C3"
color_V = "C4"

style_baseline = "-"
style_bw = "--"
style_mean = ":"

# --------------------------
# LOAD BASELINE
# --------------------------
df_base = pd.read_csv(f"./cleaned/{baseline_h_file}", delimiter=";")
power_0_h = (df_base[" S 0 [mW]"].to_numpy() + df_base[" S 1 [mW]"].to_numpy()) / 2
avg_0_h = np.mean(power_0_h)
err_0_h = np.std(power_0_h)

print(f"Baseline Power in H = {avg_0_h * 1000:.1f} ± {err_0_h * 1000:.1f} µW\n")

df_base = pd.read_csv(f"./cleaned/{baseline_v_file}", delimiter=";")
power_0_v = (df_base[" S 0 [mW]"].to_numpy() - df_base[" S 1 [mW]"].to_numpy()) / 2
avg_0_v = np.mean(power_0_v)
err_0_v = np.std(power_0_v)

print(f"Baseline Power in V = {avg_0_v * 1000:.1f} ± {err_0_v * 1000:.1f} µW\n")

files = {
    "1 BW": ("1BW_h.csv", "1BW_v.csv"),
    "2 BW": ("2BW_h.csv", "2BW_v.csv"),
    "3 BW": ("3BW_h.csv", "3BW_v.csv"),
    "4 BW": ("4BW_h.csv", "4BW_v.csv"),
    "5 BW": ("5BW_h.csv", "5BW_v.csv"),
    "6 BW": ("6BW_h.csv", "6BW_v.csv"),
    "7 BW": ("7BW_h.csv", "7BW_v.csv"),
    "8 BW": ("8BW_h.csv", "8BW_v.csv"),
}


fig, axs = plt.subplots(4, 2, figsize=(12, 10))
axs = axs.flatten()

ax_index = 0

for amount_of_bw, (file_h, file_v) in files.items():
    # --- Load H ---
    df_h = pd.read_csv(f"./cleaned/{file_h}", delimiter=";")
    ph = (df_h[" S 0 [mW]"].to_numpy() + df_h[" S 1 [mW]"].to_numpy()) / 2
    avg_h = np.mean(ph)
    err_h = np.std(ph)

    # --- Load V ---
    df_v = pd.read_csv(f"./cleaned/{file_v}", delimiter=";")
    pv = (df_v[" S 0 [mW]"].to_numpy() - df_v[" S 1 [mW]"].to_numpy()) / 2
    avg_v = np.mean(pv)
    err_v = np.std(pv)

    # --- Transmittance (H) ---
    T_h = avg_h / avg_0_h
    dT_h = T_h * np.sqrt((err_h / avg_h) ** 2 + (err_0_h / avg_0_h) ** 2)

    # --- Transmittance (V) ---
    T_v = avg_v / avg_0_v
    dT_v = T_v * np.sqrt((err_v / avg_v) ** 2 + (err_0_v / avg_0_v) ** 2)

    # --------------------------
    # PRINT RESULTS
    # --------------------------

    console = Console()

    table = Table(show_header=False, box=None)

    console.print(f"[bold cyan]--- {amount_of_bw} ---[/bold cyan]")

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
        f"{T_h * 100:.1f}% ± {dT_h * 100:.1f}%",
        f"{avg_h * 1000:.1f} ± {err_h * 1000:.1f} µW",
    )

    # V row
    t.add_row(
        "[bold red]V[/bold red]",
        f"{T_v * 100:.1f}% ± {dT_v * 100:.1f}%",
        f"{avg_v * 1000:.1f} ± {err_v * 1000:.1f} µW",
    )

    console.print(t)
    console.print("[dim]" + "─" * 42 + "[/dim]")

    ax = axs[ax_index]

    x_h = np.arange(len(ph))
    x_v = np.arange(len(pv))
    x_b_h = np.arange(len(power_0_h))
    x_b_v = np.arange(len(power_0_v))

    ax.plot(
        x_b_h,
        power_0_h,
        label="Baseline H",
        color=color_H_base,
        linestyle=style_baseline,
        linewidth=2,
        alpha=0.7,
    )

    ax.plot(
        x_b_v,
        power_0_v,
        label="Baseline V",
        color=color_V_base,
        linestyle=style_baseline,
        linewidth=2,
        alpha=0.7,
    )

    ax.plot(
        x_h,
        ph,
        label=f"{amount_of_bw} – H",
        color=color_H,
        linestyle=style_bw,
        linewidth=2,
    )

    ax.plot(
        x_v,
        pv,
        label=f"{amount_of_bw} – V",
        color=color_V,
        linestyle=style_bw,
        linewidth=2,
    )

    ax.axhline(avg_0_h, color=color_H_base, linestyle=style_mean, alpha=0.5)
    ax.axhline(avg_0_v, color=color_V_base, linestyle=style_mean, alpha=0.5)
    ax.axhline(avg_h, color=color_H, linestyle=style_mean, alpha=0.9)
    ax.axhline(avg_v, color=color_V, linestyle=style_mean, alpha=0.9)

    ax.set_title(f"{amount_of_bw}: Power vs Time")
    ax.set_xlabel("Sample index")
    ax.set_ylabel("Power [mW]")
    ax.legend(loc="center right")
    ax.grid(True)

    ax_index += 1

plt.tight_layout()
plt.show()
