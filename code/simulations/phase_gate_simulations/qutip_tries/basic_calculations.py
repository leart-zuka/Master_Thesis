from typing import Dict, Literal, Tuple, TypedDict
import numpy as np
import matplotlib.pyplot as plt
from helpers.printing import print_data, display_fidelity
from helpers.plotting import plot_two_matrices_styled_grid
from helpers.generic_computations import normalize_matrix
from rich import print

attenuate_light = True
attenuate_light = False


class params_type(TypedDict):
    mu_rf: float
    mu_fc: float
    kappa: float
    kappa_oc: float
    d_w_r: float
    d_w_a: float
    gamma: float
    g: float


def reflection_coefficient(
    mu_rf: float,
    mu_fc: float,
    mu_fc_phi: float,
    kappa: float,
    kappa_oc: float,
    d_w_r: float,
    d_w_a: float,
    gamma: float,
    g: float,
) -> complex:
    """
    Taken from Manuels's thesis (page 26)
    """
    C_c = g**2 / (2 * (gamma + 1j * d_w_a) * (kappa + 1j * d_w_r))
    r_r = mu_rf - (mu_fc * np.exp(1j * mu_fc_phi)) ** 2 * 2 * kappa_oc / (
        kappa + 1j * d_w_r
    ) * 1 / (2 * C_c + 1)
    return r_r


# def get_reflection_coefficient_error(
#     mu_rf: float,
#     mu_fc: float,
#     mu_fc_phi: float,
#     kappa: float,
#     kappa_oc: float,
#     d_w_r: float,
#     d_w_a: float,
#     gamma: float,
#     g: float,
#
#         ):
#


def phase_shift(reflection_coefficient: complex) -> np.float32:
    """
    Taken from Dominic's thesis (page 18)
    """
    phase_shift = np.angle(reflection_coefficient)
    return phase_shift


def compute_params(
    mu_rf: float,
    mu_fc: float,
    mu_fc_phi: float,
    kappa: float,
    kappa_oc: float,
    d_w_r: float,
    d_w_a: float,
    gamma: float,
    g: float,
) -> Tuple[complex, np.float32]:
    reflection_amplitude = reflection_coefficient(
        mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g
    )
    phase = phase_shift(reflection_amplitude)
    return reflection_amplitude, phase


def compute_process_fidelity(measured_matrix: np.ndarray):
    I = np.array([[1, 0], [0, 1]])
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])

    pauli_operators = [I, X, Y, Z]

    basis_of_operators = [
        np.kron(b1, b2) / 2 for b1 in pauli_operators for b2 in pauli_operators
    ]

    ideal_cnot = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])

    process_fidelity = 0
    dimensions = 4

    for operator in basis_of_operators:
        e_u = ideal_cnot @ operator @ ideal_cnot.conj().T
        e = measured_matrix @ operator @ measured_matrix.conj().T

        process_fidelity += np.trace(e_u.conj().T @ e).real
    return process_fidelity / dimensions**2


def plot_cphase_and_cnot(
    basis: Dict[
        Literal["|0,pi>", "|1,pi>", "|0,V>", "|1,V>"],
        Dict[Literal["Delta a", "Delta c"], float],
    ],
    title: str,
    attenuate_light: bool = False,
):
    results = {}

    for label, detunings in basis.items():
        params = [
            mu_rf,
            mu_fc,
            mu_fc_phi,
            kappa,
            kappa_oc,
            detunings["Delta c"],
            detunings["Delta a"],
            gamma,
            g,
        ]
        reflection_amplitude, phase = compute_params(*params)
        magnitude = np.abs(reflection_amplitude)
        power = magnitude**2
        results[label] = {
            "r": reflection_amplitude,
            "Δr": 0,
            "|r|": magnitude,
            "|r|^2": power,
            "phase_rad": phase,
            "phase_deg": np.degrees(phase),
        }

    print_data(title, basis, results)

    r_0_pi = results["|0,pi>"]["|r|^2"]
    r_1_pi = results["|1,pi>"]["|r|^2"]
    r_0_V = results["|0,V>"]["|r|^2"]
    r_1_V = results["|1,V>"]["|r|^2"]

    c_0_pi = results["|0,pi>"]["r"]
    c_1_pi = results["|1,pi>"]["r"]
    c_0_V = results["|0,V>"]["r"]
    c_1_V = results["|1,V>"]["r"]

    if attenuate_light:
        avg_pi = (r_0_pi + r_1_pi) / 22
        avg_v = (r_0_V + r_1_V) / 22
        reduction_v_polarization = avg_pi / avg_v

        r_0_V *= reduction_v_polarization
        r_1_V *= reduction_v_polarization
        c_0_V *= np.sqrt(reduction_v_polarization)
        c_1_V *= np.sqrt(reduction_v_polarization)

    input_matrix = np.array(
        [
            [r_0_pi, 0.0, 0.0, 0.0],
            [0.0, r_1_pi, 0.0, 0.0],
            [0.0, 0.0, r_0_V, 0.0],
            [0.0, 0.0, 0.0, r_1_V],
        ]
    )
    cnot_matrix = np.array(
        [
            [
                np.abs((c_1_pi + c_1_V) / 2) ** 2,
                np.abs((-c_1_pi + c_1_V) / 2) ** 2,
                0,
                0,
            ],
            [
                np.abs((-c_1_pi + c_1_V) / 2) ** 2,
                np.abs((c_1_pi + c_1_V) / 2) ** 2,
                0,
                0,
            ],
            [
                0,
                0,
                np.abs((+c_0_pi + c_0_V) / 2) ** 2,
                np.abs((-c_0_pi + c_0_V) / 2) ** 2,
            ],
            [
                0,
                0,
                np.abs((-c_0_pi + c_0_V) / 2) ** 2,
                np.abs((+c_0_pi + c_0_V) / 2) ** 2,
            ],
        ]
    )

    normalized_input = normalize_matrix(input_matrix)
    normalized_cnot = normalize_matrix(cnot_matrix)

    c_phase_labels = ["|0,π⟩", "|1,π⟩", "|0,V⟩", "|1,V⟩"]
    cnot_labels = ["|1,+⟩", "|1,-⟩", "|0,+⟩", "|0,-⟩"]

    plot_two_matrices_styled_grid(
        mat_left=input_matrix,
        norm_left=normalized_input,
        labels_left=c_phase_labels,
        title_left="Reflection Intensity (C-PHASE)",
        mat_right=cnot_matrix,
        norm_right=normalized_cnot,
        labels_right=cnot_labels,
        title_right="Population Probability (CNOT)",
        title_figure=title,
    )
    return normalized_cnot


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
g = 2 * np.pi * 0.0386  # 2-1' splitting
kappa = 2 * np.pi * 0.058
kappa_oc = kappa * 0.85
gamma = 2 * np.pi * 0.006065
mu_rf = 0.978
mu_fc = 0.873
mu_fc_phi = 0.024

# Errors
g_sigma = 2 * np.pi * 0.004
kappa_sigma = 2 * np.pi * 0.00037 / 2
kappa_oc_sigma = kappa_sigma * 0.085
gamma_sigma = 2 * np.pi * 0.000018
mu_rf_sigma = 0.006
mu_fc_sigma = 0.002
mu_fc_phi_sigma = 0.0001

# Reduction for polarization
reduction_v_pol = 0.362

if attenuate_light:
    theoretical_cnot = plot_cphase_and_cnot(basis, "2-1' Transition", reduction_v_pol)
else:
    theoretical_cnot = plot_cphase_and_cnot(basis, "2-1' Transition")

process_fidelity = compute_process_fidelity(theoretical_cnot)
display_fidelity(process_fidelity)
#
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
