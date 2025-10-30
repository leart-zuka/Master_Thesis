from typing import Dict, Literal
import numpy as np
from helpers.printing import display_fidelity
from helpers.plotting import plot_cphase_and_cnot
from helpers.compute_reflection_parameters import (
    compute_process_fidelity,
    params_type,
)

attenuate_light = True
# attenuate_light = False

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
    "g": 2 * np.pi * 0.0353,
    "kappa": 2 * np.pi * 0.058,
    "kappa_oc": 2 * np.pi * 0.058 * 0.85,
    "gamma": 2 * np.pi * 0.006065,
    "mu_rf": 0.978,
    "mu_fc": 0.873,
    "mu_fc_phi": 0.024,
}

# Errors
params_err_dir: params_type = {
    "g": 2 * np.pi * 0.0036,
    "kappa": 2 * np.pi * 0.00037 / 2,
    "kappa_oc": 2 * np.pi * 0.00037 / 2 * 0.85,
    "gamma": 2 * np.pi * 0.000018,
    "mu_rf": 0.006,
    "mu_fc": 0.002,
    "mu_fc_phi": 0.001,
}

# Reduction for polarization
reduction_v_pol = 0.362


if attenuate_light:
    theoretical_c_phase, theoretical_cnot = plot_cphase_and_cnot(
        basis, "2-1' Transition", params_dir, params_err_dir, attenuate_light
    )
else:
    theoretical_c_phase, theoretical_cnot = plot_cphase_and_cnot(
        basis, "2-1' Transition", params_dir, params_err_dir
    )

process_fidelity = compute_process_fidelity(theoretical_cnot)
display_fidelity(process_fidelity)

# g = 2 * np.pi * 0.0662  # 2-2' splitting
# Kappa = 2 * np.pi * 0.058
# Kappa_oc = Kappa * 0.85
# Gamma = 2 * np.pi * 0.006065 / 2
# Coopoerativity = g**2 / (2 * Kappa * Gamma)
#
# basis: Dict[
#     Literal["|0,pi>", "|1,pi>", "|0,V>", "|1,V>"],
#     Dict[Literal["Delta a", "Delta c"], float],
# ] = {
#     "|0,pi>": {
#         "Delta a": 2 * np.pi * 6.385,
#         "Delta c": 2 * np.pi * 0,
#     },
#     "|1,pi>": {
#         "Delta a": 2 * np.pi * 0,
#         "Delta c": 2 * np.pi * 0,
#     },
#     "|0,V>": {
#         "Delta a": 2 * np.pi * 6.385,
#         "Delta c": 2 * np.pi * 0.5,
#     },
#     "|1,V>": {
#         "Delta a": 2 * np.pi * 0,
#         "Delta c": 2 * np.pi * 0.5,
#     },
# }
#
# reduction_v_pol = 0.465
#
# if attenuate_light:
#     theoretical_cnot = plot_cphase_and_cnot(basis, "2-2' Transition", reduction_v_pol)
# else:
#     theoretical_cnot = plot_cphase_and_cnot(basis, "2-2' Transition")
#
#
# process_fidelity = compute_process_fidelity(theoretical_cnot)
# display_fidelity(process_fidelity)
