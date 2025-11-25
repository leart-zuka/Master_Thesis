import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8")

# --------------------------
# FILES
# --------------------------
files = {
    "BW1": ("bw1_h.csv", "bw1_v.csv"),
    "BW2": ("bw2_h.csv", "bw2_v.csv"),
}

baseline_file = "baseline.csv"

# --------------------------
# LOAD BASELINE
# --------------------------
df_base = pd.read_csv(f"./cleaned/{baseline_file}", delimiter=";")
power_0 = df_base[" S 0 [mW]"].to_numpy()
avg_0 = np.mean(power_0)
err_0 = np.std(power_0) / np.sqrt(len(power_0))

print(f"Baseline Power = {avg_0 * 1000:.3f} ± {err_0 * 1000:.3f} µW\n")


# --------------------------
# PLOTTING GRID
# --------------------------
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

# Axis index
ax_index = 0


# --------------------------
# PROCESS EACH BREWSTER WINDOW
# --------------------------
for bw_name, (file_h, file_v) in files.items():
    # --- Load H ---
    df_h = pd.read_csv(f"./cleaned/{file_h}", delimiter=";")
    ph = df_h[" S 0 [mW]"].to_numpy()
    avg_h = np.mean(ph)
    err_h = np.std(ph) / np.sqrt(len(ph))

    # --- Load V ---
    df_v = pd.read_csv(f"./cleaned/{file_v}", delimiter=";")
    pv = df_v[" S 0 [mW]"].to_numpy()
    avg_v = np.mean(pv)
    err_v = np.std(pv) / np.sqrt(len(pv))

    # --- Transmittance (H) ---
    T_h = avg_h / avg_0
    dT_h = T_h * np.sqrt((err_h / avg_h) ** 2 + (err_0 / avg_0) ** 2)

    # --- Transmittance (V) ---
    T_v = avg_v / avg_0
    dT_v = T_v * np.sqrt((err_v / avg_v) ** 2 + (err_0 / avg_0) ** 2)

    # --------------------------
    # PRINT RESULTS
    # --------------------------
    print(f"--- {bw_name} ---")
    print(f"H: {T_h * 100:.2f}% ± {dT_h * 100:.2f}%")
    print(f"V: {T_v * 100:.2f}% ± {dT_v * 100:.2f}%")
    print("------------------------------------------")

    # --------------------------
    # PLOT
    # --------------------------
    ax = axs[ax_index]

    # x-axis indices just to visualize raw power time series
    x_h = np.arange(len(ph))
    x_v = np.arange(len(pv))
    x_b = np.arange(len(power_0))

    ax.plot(x_b, power_0, label="Baseline", linewidth=2)
    ax.plot(x_h, ph, label=f"{bw_name} – H")
    ax.plot(x_v, pv, label=f"{bw_name} – V")

    # Mark average lines
    ax.axhline(avg_0, color="black", linestyle="--", alpha=0.6)
    ax.axhline(avg_h, color="C1", linestyle="--", alpha=0.6)
    ax.axhline(avg_v, color="C2", linestyle="--", alpha=0.6)

    ax.set_title(f"{bw_name}: Power vs Time")
    ax.set_xlabel("Sample index")
    ax.set_ylabel("Power [mW]")
    ax.legend(loc="upper right")
    ax.grid(True)

    ax_index += 1


plt.tight_layout()
plt.show()
