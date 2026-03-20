from typing import Dict, Literal
import numpy as np
from helpers.printing import display_fidelity
from helpers.plotting import plot_cphase_and_cnot
from helpers.compute_reflection_parameters import (
    compute_process_fidelity,
    compute_fidelity_from_prob_matrix,
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
    attenuate_light,
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

    process_fidelity = compute_fidelity_from_prob_matrix(
        gate_matrices["normalized_cnot"], basis="cnot"
    )
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

plt.figure(figsize=(12, 8))
plt.plot(x, y, linewidth=2)
plt.xlabel("Total Transmission for V polarized light", fontsize=24)
plt.ylabel(r"$F_{process}$", fontsize=24)
plt.tick_params(axis="both", which="major", labelsize=20, length=10, width=2)

# Highlight max slope point
plt.scatter(x_max, y_max, s=120, color="red", zorder=5, marker="*", label="Optimal")

# Annotate with arrow + values
# Use 'offset points' for xytext to make it relative and stable
plt.annotate(
    f"T = {x_max:.3f}\n" + r"$F_{proc}$ = " + f"{y_max:.3f}",
    xy=(x_max, y_max),
    xytext=(-80, -50),
    textcoords="offset points",
    arrowprops=dict(
        arrowstyle="->", connectionstyle="arc3,rad=.2", lw=1.5, color="black"
    ),
    fontsize=10,
    fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9, ec="gray"),
)

for i, xn in enumerate(attenuations_additions):
    yn = float(np.interp(xn, x, y))
    plt.axvline(x=xn, color="coral", linestyle="--", linewidth=1, alpha=0.5)
    plt.scatter([xn], [yn], color="coral", s=40, zorder=4, alpha=0.8)

    # --- STAGGERING LOGIC ---
    # Alternate labels above and below the point to prevent overlapping at high BW counts
    direction = 1 if i % 2 == 0 else -1
    if i == 7:
        direction = 1
    y_offset = 12 * direction  # pixels
    if i == 7:
        plt.annotate(
            f"{i + 1} BWs",
            xy=(xn, yn),
            xytext=(-10, y_offset),
            textcoords="offset points",
            fontsize=14,
            fontweight="bold",
            va="bottom" if direction > 0 else "top",
        )
    else:
        plt.annotate(
            f"{i + 1} BWs",
            xy=(xn, yn),
            xytext=(5, y_offset),
            textcoords="offset points",
            fontsize=14,
            fontweight="bold",
            va="bottom" if direction > 0 else "top",
        )

# --- THE FIXES ---
# plt.title("Process Fidelity vs V-Polarized Transmission", pad=25, fontsize=12)
# plt.grid(True, alpha=0.3)

# 1. Add internal margins: 0.15 on Y ensures the staggered labels fit inside the frame
plt.margins(x=0.1, y=0.2)

# 2. Use a tighter Y limit if needed to focus on the high-fidelity region
plt.ylim(bottom=min(y) * 0.95)
plt.ylim([0.5, 1])

# 3. Final layout and save
plt.tight_layout()
plt.savefig("./plots/attenuation_fidelity.svg", bbox_inches="tight")
plt.show()
