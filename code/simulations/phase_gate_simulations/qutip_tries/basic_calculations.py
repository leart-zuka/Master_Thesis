from typing import Dict, Literal
import numpy as np
from helpers.printing import display_fidelity
from helpers.plotting import plot_cphase_and_cnot
from helpers.compute_reflection_parameters import (
    compute_process_fidelity,
    params_type,
)

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

gate_matrices = plot_cphase_and_cnot(
    basis,
    "2-1' Transition",
    params_dir,
    params_err_dir,
    attenuate_light,
    save_fig=True,
)

process_fidelity = compute_process_fidelity(gate_matrices["normalized_cnot"])
process_fidelity_err = compute_process_fidelity(gate_matrices["normalized_cnot"])
display_fidelity(process_fidelity, error=process_fidelity_err)


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
