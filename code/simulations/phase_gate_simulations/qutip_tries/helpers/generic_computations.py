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
    out_plus,
    out_minus,
    dt,
    p_dc=1e-4,
    eta=0.74205,
    pulse_length: float = 1e-6,
    detector_dead_time: float = 20e-9,
    test: bool = False,
):

    n_out_plus = np.sum(np.abs(out_plus) ** 2) * dt
    n_out_minus = np.sum(np.abs(out_minus) ** 2) * dt

    n_tot = n_out_plus + n_out_minus

    l_plus = n_out_plus * eta
    l_minus = n_out_minus * eta

    gamma_p = l_plus + p_dc
    gamma_m = l_minus + p_dc

    # Signal-to-noise ratio component
    P_SNR = (l_plus + l_minus) / (l_plus + l_minus + 2 * p_dc)

    P_SNR_p = (l_plus * np.exp(-(n_out_plus + p_dc))) / (1 - np.exp(-(l_plus + p_dc)))
    P_SNR_m = (l_minus * np.exp(-(n_out_minus + p_dc))) / (
        1 - np.exp(-(l_minus + p_dc))
    )

    prob_distinction = (
        2 * pulse_length * detector_dead_time - pulse_length**2
    ) / pulse_length**2

    P_click_p = (
        gamma_p * np.exp(-gamma_p)
        + (gamma_p) ** 2 / 2 * np.exp(-gamma_p) * prob_distinction
    )
    P_click_m = (
        gamma_m * np.exp(-gamma_m)
        + (gamma_m) ** 2 / 2 * np.exp(-gamma_m) * prob_distinction
    )

    P_SNR_p = l_plus * np.exp(-gamma_p) / P_click_p
    P_SNR_m = l_minus * np.exp(-gamma_m) / P_click_m

    R_plus = n_out_plus / (n_out_plus + n_out_minus)
    R_minus = n_out_minus / (n_out_plus + n_out_minus)

    # Calculate final operator probabilities
    if not test:
        P_out_plus = (P_SNR_p * R_plus) + ((1 - P_SNR_p) * 0.5)
        P_out_minus = (P_SNR_m * R_minus) + ((1 - P_SNR_m) * 0.5)
    else:
        P_out_plus = (P_SNR_p * R_plus) + ((1 - P_SNR_p) * 0.5)
        P_out_minus = (P_SNR_p * R_minus) + ((1 - P_SNR_m) * 0.5)

    return P_out_plus, P_out_minus
