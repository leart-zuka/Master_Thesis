import matplotlib.pyplot as plt
import numpy as np

n_bw = np.arange(1, 9)
T_H = np.array([99.96, 99.94, 99.15, 99.20, 99.05, 99.06, 98.90, 98.86])
T_H_err = np.array([0.18, 0.17, 0.17, 0.18, 0.17, 0.18, 0.16, 0.16])
T_V = np.array([77.38, 61.10, 52.29, 45.09, 31.82, 26.80, 22.47, 16.56])
T_V_err = np.array([0.13, 0.12, 0.10, 0.09, 0.10, 0.06, 0.04, 0.04])

fig, ax = plt.subplots(figsize=(5.5, 3.8))

ax.errorbar(
    n_bw,
    T_H,
    yerr=T_H_err,
    fmt="o-",
    color="#2166ac",
    capsize=3,
    markersize=5,
    linewidth=1.2,
    label=r"$T_H$",
)
ax.errorbar(
    n_bw,
    T_V,
    yerr=T_V_err,
    fmt="s-",
    color="#b2182b",
    capsize=3,
    markersize=5,
    linewidth=1.2,
    label=r"$T_V$",
)

target_TV = 24.0
ax.axhline(
    target_TV,
    color="green",
    linestyle="--",
    linewidth=0.8,
    label=f"Target $T_V = {target_TV}\%$",
)

ax.set_xlabel("Number of Brewster windows", fontsize=13)
ax.set_ylabel("Transmission (%)", fontsize=13)
ax.set_xticks(n_bw)
ax.set_ylim(0, 105)
ax.legend(frameon=False, fontsize=12)
ax.tick_params(
    axis="both", labelsize=12, direction="in", which="both", top=True, right=True
)

fig.tight_layout()
fig.savefig("bw_transmission_presentation.pdf", bbox_inches="tight")
fig.savefig("bw_transmission_presentation.svg", bbox_inches="tight")
fig.savefig("bw_transmission_presentation.png", dpi=200, bbox_inches="tight")
print("Plot saved.")
