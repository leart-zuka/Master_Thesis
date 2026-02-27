from typing import Dict, Literal
import numpy as np
from helpers.printing import display_fidelity
from helpers.plotting import plot_cphase_and_cnot
from helpers.compute_reflection_parameters import (
    compute_process_fidelity,
    params_type,
)
import matplotlib.pyplot as plt

attenuate_light = True

basis: Dict[
    Literal["|0,pi>", "|1,pi>", "|0,V>", "|1,V>"],
    Dict[Literal["Delta a", "Delta c"], float],
] = {
    "|0,pi>": {
        "Delta a": 2 * np.pi * 6.385,
        "Delta c": 2 * np.pi * 0,
    },
    "|1,pi>": {
        "Delta a": 2 * np.pi * 0,
        "Delta c": 2 * np.pi * 0,
    },
    "|0,V>": {
        "Delta a": 2 * np.pi * 6.385,
        "Delta c": 2 * np.pi * 0.5,
    },
    "|1,V>": {
        "Delta a": 2 * np.pi * 0,
        "Delta c": 2 * np.pi * 0.5,
    },
}

# Values
params_dir: params_type = {
    "g": 2 * np.pi * 0.024,
    "kappa": 2 * np.pi * 0.058,
    "kappa_oc": 2 * np.pi * 0.058 * 0.85,
    "gamma": 2 * np.pi * 0.006065 / 2,
    "mu_rf": 0.978,
    "mu_fc": 0.873,
    "mu_fc_phi": 0.024,
}

# Errors
params_err_dir: params_type = {
    "g": 2 * np.pi * 0.0007,
    "kappa": 2 * np.pi * 0.00037 / 2,
    "kappa_oc": 2 * np.pi * 0.00037 / 2 * 0.85,
    "gamma": 2 * np.pi * 0.000018,
    "mu_rf": 0.006,
    "mu_fc": 0.002,
    "mu_fc_phi": 0.001,
}

gate_matrices = plot_cphase_and_cnot(
    basis,
    "2-1' Transition",
    params_dir,
    params_err_dir,
    # attenuate_light,
    special_attenuation=0.24,
    save_fig=True,
)

process_fidelity = compute_process_fidelity(gate_matrices["normalized_cnot"])
process_fidelity_err = compute_process_fidelity(gate_matrices["cnot_matrix_err"])
display_fidelity(process_fidelity, error=process_fidelity_err)


transmissions = np.linspace(0.01, 1, 100)
fidelities = []

for transmission in transmissions:
    gate_matrices = plot_cphase_and_cnot(
        basis,
        "2-1' Transition",
        params_dir,
        params_err_dir,
        attenuate_light,
        special_attenuation=transmission,
        display_data=False,
        show_fig=False,
        save_fig=False,
    )

    process_fidelity = compute_process_fidelity(gate_matrices["normalized_cnot"])
    process_fidelity_err = compute_process_fidelity(gate_matrices["normalized_cnot"])
    fidelities.append(process_fidelity)


x = transmissions
y = np.array(fidelities)

attenuations = [0.7738, 0.7736, 0.7754, 0.7746, 0.7735, 0.7716, 0.7740, 0.7753]

attenuations_additions = np.cumprod(attenuations).tolist()
attenuations_additions = [
    0.7738,
    0.6110,
    0.5229,
    0.4509,
    0.3182,
    0.2680,
    0.2247,
    0.1656,
]
# Find max index
idx = np.argmax(y)
x_max = x[idx]
y_max = y[idx]

# Plot
# Plot
plt.figure(figsize=(7, 5))  # Increased height slightly to accommodate title + labels
plt.plot(x, y)
plt.xlabel("Total Transmission for V polarized light")
plt.ylabel(r"$F_{process}$")

# Highlight max slope point
plt.scatter(x_max, y_max, s=80, color="red", zorder=5, marker="*")

# Annotate with arrow + values
# Improved xytext logic to prevent clipping on the left
plt.annotate(
    f"T = {x_max:.3f}\n" + r"$F_{proc}$ = " + f"{y_max:.3f}",
    xy=(x_max, y_max),
    xytext=(x_max - 0.15, y_max - 0.1),
    arrowprops=dict(arrowstyle="->", lw=1.5, color="black"),
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8, ec="gray"),
)

for i, xn in enumerate(attenuations_additions):
    yn = float(np.interp(xn, x, y))
    plt.axvline(x=xn, color="coral", linestyle="--", linewidth=1, alpha=0.7)
    plt.scatter([xn], [yn], color="coral", s=40, zorder=4, alpha=0.7)

    plt.annotate(
        f"{i + 1} BWs",
        xy=(xn, yn),
        xytext=(5, 5),  # 5pt offset from the point
        textcoords="offset points",
        fontsize=9,
        fontweight="bold",
    )

# --- THE FIXES ---
plt.title("Process Fidelity vs V-Polarized Transmission", pad=20)
plt.grid(True, alpha=0.3)

# 1. Add internal margins so labels near the edges aren't cut off
plt.margins(x=0.1, y=0.15)

# 2. Use constrained_layout or one tight_layout call
plt.tight_layout()

# 3. Use bbox_inches="tight" during save to capture everything
plt.savefig("./plots/attenuation_fidelity.svg", bbox_inches="tight")
plt.show()
