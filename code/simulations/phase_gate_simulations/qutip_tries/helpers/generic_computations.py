import numpy as np


def normalize_matrix(input_matrix: np.ndarray):
    """
    note: make sure that all entries are in terms of probabilities and not just raw complex numbers
          as that can lead to the normalization costant being close to 0, ergo blowing up your term
          which you're trying to normalize :>
    """
    normalized = np.zeros_like(input_matrix)
    normalization_constants = normalizations(input_matrix)
    normalized = input_matrix / normalization_constants
    return normalized


def normalizations(input_matrix: np.ndarray):
    col_sums = input_matrix.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    return np.abs(col_sums)


def calculate_detection_probabilities(
    n_out_plus,
    n_out_minus,
    dt,
    r_dc=1e-4,
    eta=0.74205,
    pulse_length: float = 1e-6,
    detector_dead_time: float = 20e-9,
):

    n_out_p = np.sum(n_out_plus) * dt
    n_out_m = np.sum(n_out_minus) * dt

    l_p = n_out_p * eta
    l_m = n_out_m * eta

    gamma_p = l_p + r_dc
    gamma_m = l_m + r_dc

    # Probability that two events happen within the detector dead time
    prob_distinction = (
        2 * pulse_length * detector_dead_time - detector_dead_time**2
    ) / pulse_length**2

    # Probability that exactly 1 click was produced, consists of
    # two cases:
    # actually one event
    # + two events that happen within detector dead-time
    F_photons = (l_p + l_m) / (
        gamma_p * (1 + (gamma_p / 2) * prob_distinction)
        + gamma_m * (1 + (gamma_m / 2) * prob_distinction)
    )

    R_p = n_out_p / (n_out_p + n_out_m)
    R_m = n_out_m / (n_out_p + n_out_m)

    # Calculate final operator probabilities
    P_out_plus = F_photons * R_p + (1 - F_photons) * 0.5
    P_out_minus = F_photons * R_m + (1 - F_photons) * 0.5

    return P_out_plus, P_out_minus
