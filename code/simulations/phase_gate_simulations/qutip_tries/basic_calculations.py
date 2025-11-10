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
    "gamma": 2 * np.pi * 0.006065,
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
    special_attenuation=0.17,
    save_fig=True,
)

process_fidelity = compute_process_fidelity(gate_matrices["normalized_cnot"])
process_fidelity_err = compute_process_fidelity(gate_matrices["cnot_matrix_err"])
display_fidelity(process_fidelity, error=process_fidelity_err)


# transmissions = np.linspace(0.01, 1, 100)
# fidelities = []
#
# for transmission in transmissions:
#     gate_matrices = plot_cphase_and_cnot(
#         basis,
#         "2-1' Transition",
#         params_dir,
#         params_err_dir,
#         attenuate_light,
#         special_attenuation=transmission,
#         show_fig=False,
#         save_fig=False,
#     )
#
#     process_fidelity = compute_process_fidelity(gate_matrices["normalized_cnot"])
#     process_fidelity_err = compute_process_fidelity(gate_matrices["normalized_cnot"])
#     fidelities.append(process_fidelity)
#
#
# x = transmissions
# y = np.array(fidelities)
#
# attenuations = [0.7773, 0.7778, 0.7744, 0.7731]
# attenuations_additions = [
#     attenuations[0],
#     attenuations[0] * attenuations[1],
#     attenuations[0] * attenuations[1] * attenuations[2],
#     attenuations[0] * attenuations[1] * attenuations[2] * attenuations[3],
# ]
# # attenuations_additions = [
# #     1 - sum(attenuations[: i + 1]) for i in range(len(attenuations))
# # ]
#
# # Find max index
# idx = np.argmax(y)
# x_max = x[idx]
# y_max = y[idx]
#
# # Plot
# plt.figure(figsize=(7, 4))
# plt.plot(x, y)
# plt.xlabel("Total Tranmission for V polarized light")
# plt.ylabel(r"$F_{process}$")
#
# # Highlight max slope point
# plt.scatter(x_max, y_max, s=80, color="red", zorder=5)
#
# # Annotate with arrow + values
# plt.annotate(
#     f"T = {x_max:.3f}\nF_proc = {y_max:.3f}",
#     xy=(x_max, y_max),
#     xytext=(x_max - 0.3 * (max(y) - min(x)), y_max - 0.08 * (max(y) - min(y))),
#     arrowprops=dict(arrowstyle="->", lw=1.5),
#     fontsize=10,
#     bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
# )
# #
# for xn in attenuations_additions:
#     print("attenuation (fraction):", xn)
#     yn = float(np.interp(xn, x, y))
#     print("interpolated y:", yn)
#
#     plt.axvline(x=xn, color="red", linestyle="--", linewidth=1)
#     plt.scatter([xn], [yn], color="red", s=40)
#
#     plt.annotate(
#         f"({xn:.3f}, {yn:.3f})",
#         xy=(xn, yn),
#         xytext=(xn + 0.005, yn + 0.01),
#         fontsize=9,
#     )
#
# plt.grid(True, alpha=0.3)
# plt.tight_layout()
# plt.title("Process Fidelity for different Transmission of V polarized light")
# plt.tight_layout()
# plt.savefig("attenuation_fidelity.svg")
# plt.show()
#
params_dir: params_type = {
    "g": 2 * np.pi * 0.0642,
    "kappa": 2 * np.pi * 0.058,
    "kappa_oc": 2 * np.pi * 0.058 * 0.85,
    "gamma": 2 * np.pi * 0.006065,
    "mu_rf": 0.978,
    "mu_fc": 0.873,
    "mu_fc_phi": 0.024,
}

# Errors
params_err_dir: params_type = {
    "g": 2 * np.pi * 0.0026,
    "kappa": 2 * np.pi * 0.00037 / 2,
    "kappa_oc": 2 * np.pi * 0.00037 / 2 * 0.85,
    "gamma": 2 * np.pi * 0.000018,
    "mu_rf": 0.006,
    "mu_fc": 0.002,
    "mu_fc_phi": 0.001,
}

# theoretical_c_phase, theoretical_cnot, theoretical_cnot_err = plot_cphase_and_cnot(
#     basis,
#     "2-2' Transition",
#     params_dir,
#     params_err_dir,
#     attenuate_light,
#     save_fig=True,
# )

# process_fidelity = compute_process_fidelity(theoretical_cnot)
# process_fidelity_err = compute_process_fidelity(theoretical_cnot_err)
# display_fidelity(process_fidelity, process_fidelity_err)
