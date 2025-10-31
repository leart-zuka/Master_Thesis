import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

files = [
    "bw1_h_lipo.csv",
    "bw1_v_lipo.csv",
    "bw2_h_lipo.csv",
    "bw2_v_lipo.csv",
    "bw3_h_lipo.csv",
    "bw3_v_lipo.csv",
    "bw4_h_lipo.csv",
    "bw4_v_lipo.csv",
]

for i in range(0, len(files) - 1, 2):
    df_h = pd.read_csv(files[i], delimiter=";")
    power_h = df_h[" S 0 [mW]"].to_numpy()
    average_power_h = np.average(power_h)
    std_power_h = np.std(power_h) / np.sqrt(len(power_h))

    print(
        f"For H polarization: {average_power_h * 1000:.3f} +/- {std_power_h * 1000:.3f} μW"
    )

    df_v = pd.read_csv(files[i + 1], delimiter=";")
    power_v = df_v[" S 0 [mW]"].to_numpy()
    average_power_v = np.average(power_v)
    std_power_v = np.std(power_v) / np.sqrt(len(power_v))

    print(
        f"For V polarization: {average_power_v * 1000:.3f} +/- {std_power_v * 1000:.3f} μW"
    )

    # Ratio + error propagation
    reduction_ratio = average_power_v / average_power_h
    std_reduction_ratio = reduction_ratio * np.sqrt(
        (std_power_v / average_power_v) ** 2 + (std_power_h / average_power_h) ** 2
    )

    # Convert to percentage (reduction = 1 - ratio)
    reduction_percent = (1 - reduction_ratio) * 100
    std_reduction_percent = std_reduction_ratio * 100

    print(
        f"Reduction for V polarization from Brewster Window {files[i][2]}: "
        f"{reduction_percent:.2f}% +/- {std_reduction_percent:.2f}%"
    )
