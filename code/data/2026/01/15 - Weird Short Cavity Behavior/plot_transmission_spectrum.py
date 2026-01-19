import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv(
    "./C1Short_Cavity_Intracavity_Trap00000.csv",
    sep=",",
    skiprows=4,
)

plt.figure(figsize=(9, 5))

plt.plot(
    df["Time"],
    df["Ampl"],
    linewidth=2.0,
    alpha=0.9,
)

plt.xlabel("Time", fontsize=13)
plt.ylabel("Amplitude", fontsize=13)
plt.title("Intracavity Trap Signal Scan", fontsize=15)

plt.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.6)

plt.tick_params(axis="both", labelsize=11)

plt.tight_layout()
plt.show()
