from typing import Dict, Tuple
import numpy as np
from helpers.printing import print_data
from helpers.plotting import plot_raw_and_normalized_styled
from helpers.generic_computations import normalize_matrix
import matplotlib.pyplot as plt
from rich import print


def reflection_coefficient(
    mu_rf: float,
    mu_fc: float,
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
    r_r = mu_rf - mu_fc**2 * 2 * kappa_oc / ((kappa + 1j * d_w_r) * (2 * C_c + 1))
    return r_r


def phase_shift(reflection_coefficient: complex) -> np.float32:
    """
    Taken from Dominic's thesis (page 18)
    """
    phase_shift = np.angle(reflection_coefficient)
    return phase_shift


def compute_params(
    mu_rf: float,
    mu_fc: float,
    kappa: float,
    kappa_oc: float,
    d_w_r: float,
    d_w_a: float,
    gamma: float,
    g: float,
) -> Tuple[complex, np.float32]:
    reflection_amplitude = reflection_coefficient(
        mu_rf, mu_fc, kappa, kappa_oc, d_w_r, d_w_a, gamma, g
    )
    phase = phase_shift(reflection_amplitude)
    return reflection_amplitude, phase


G0_kc = 2 * np.pi * 0.0386  # 2-1' splitting
# G0_kc *= np.sqrt(30 / 24)
Kappa = 2 * np.pi * 0.058
Kappa_oc = Kappa * 0.85
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2
Coopoerativity = G0_kc**2 / (2 * Kappa * Gamma_5P32_5S)
Mu_rf = 0.88
Mu_fc = 0.88

# --------------------
basis: Dict[str, Dict[str, float]] = {
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

results = {}

for label, detunings in basis.items():
    reflection_amplitude, phase = compute_params(
        mu_rf=Mu_rf,
        mu_fc=Mu_fc,
        kappa=Kappa,
        kappa_oc=Kappa_oc,
        d_w_r=basis[label]["Delta c"],
        d_w_a=basis[label]["Delta a"],
        gamma=Gamma_5P32_5S,
        g=G0_kc,
    )
    magnitude = np.abs(reflection_amplitude)
    power = magnitude**2
    results[label] = {
        "r": reflection_amplitude,
        "|r|": magnitude,
        "|r|^2": power,
        "phase_rad": phase,
        "phase_deg": np.degrees(phase),
    }

print_data(basis, results)

c_0_pi = results["|0,pi>"]["|r|^2"]
c_1_pi = results["|1,pi>"]["|r|^2"]
c_0_V = results["|0,V>"]["|r|^2"]
c_1_V = results["|1,V>"]["|r|^2"]

c_0_V *= 0.46
c_1_V *= 0.46

input_matrix = np.array(
    [
        [c_0_pi, 0.0, 0.0, 0.0],
        [0.0, c_1_pi, 0.0, 0.0],
        [0.0, 0.0, c_0_V, 0.0],
        [0.0, 0.0, 0.0, c_1_V],
    ],
    dtype=np.float32,
)
normalized_input = normalize_matrix(input_matrix)

labels = ["|0,π⟩", "|1,π⟩", "|0,V⟩", "|1,V⟩"]
plot_raw_and_normalized_styled(
    input_matrix,
    normalized_matrix=normalized_input,
    labels=labels,
    title="Reflection Intensity",
)


labels = ["|1,+⟩", "|1,-⟩", "|0,+⟩", "|0,-⟩"]

cnot_matrix = np.array(
    [
        [np.abs((c_1_pi + c_1_V) / 2) ** 2, np.abs((-c_1_pi + c_1_V) / 2) ** 2, 0, 0],
        [np.abs((-c_1_pi + c_1_V) / 2) ** 2, np.abs((c_1_pi + c_1_V) / 2) ** 2, 0, 0],
        [0, 0, np.abs((-c_0_pi + c_0_V) / 2) ** 2, np.abs((c_0_pi + c_0_V) / 2) ** 2],
        [0, 0, np.abs((c_0_pi + c_0_V) / 2) ** 2, np.abs((-c_0_pi + c_0_V) / 2) ** 2],
    ]
)

plot_raw_and_normalized_styled(
    cnot_matrix, labels=labels, title="Population Probability"
)
