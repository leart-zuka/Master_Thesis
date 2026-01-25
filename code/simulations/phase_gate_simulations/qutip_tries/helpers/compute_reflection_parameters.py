import numpy as np
from typing import Tuple, TypedDict

CNOT_IDEAL = np.array(
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex
)


class params_type(TypedDict):
    mu_rf: float
    mu_fc: float
    mu_fc_phi: float
    kappa: float
    kappa_oc: float
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

    ideal_cnot = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])

    process_fidelity = 0
    dimensions = 4

    for operator in basis_of_operators:
        e_u = ideal_cnot @ operator @ ideal_cnot.conj().T
        e = measured_matrix @ operator @ measured_matrix.conj().T

        process_fidelity += np.trace(e_u.conj().T @ e).real
    return process_fidelity / dimensions**2


def compute_signal_fidelity(measured_matrix: np.ndarray):
    ideal_cnot = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex
    )
    V_cols = [measured_matrix[:, j] for j in range(4)]
    U_cols = [ideal_cnot[:, j] for j in range(4)]

    F_sig_cols = np.array(
        [np.abs(np.vdot(U_cols[j], V_cols[j])) ** 2 for j in range(4)]
    )

    F_sig = 1 / 4 * np.sum(F_sig_cols)
    return F_sig


def normalize(vec):
    n = np.linalg.norm(vec)
    if n == 0:
        return vec
    return vec / n


def conditional_output_state(CNOT, input_index):
    """
    Returns normalized conditional output state |ψ_out⟩
    """
    out_state = CNOT[:, input_index]
    return normalize(out_state.astype(complex))


def ideal_output_state(input_index):
    basis = np.eye(4, dtype=complex)
    return CNOT_IDEAL @ basis[:, input_index]


def state_fidelity(psi_ideal, psi_actual):
    return np.abs(np.vdot(psi_ideal, psi_actual)) ** 2


def conditional_process_fidelity(CNOT):
    """
    Average conditional fidelity over computational basis inputs
    """
    fidelities = []

    for i in range(4):
        psi_out = conditional_output_state(CNOT, i)
        psi_ideal = ideal_output_state(i)

        F = state_fidelity(psi_ideal, psi_out)
        fidelities.append(F)

    return np.mean(fidelities)
