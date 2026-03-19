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
    out_plus, out_minus, dt, p_dc=1e-4, eta=0.74205, test: bool = False
):

    n_out_plus = np.sum(np.abs(out_plus) ** 2) * dt
    n_out_minus = np.sum(np.abs(out_minus) ** 2) * dt

    n_tot = n_out_plus + n_out_minus

    l_plus = n_out_plus * eta
    l_minus = n_out_minus * eta

    # Signal-to-noise ratio component
    P_SNR = (l_plus + l_minus) / (l_plus + l_minus + 2 * p_dc)

    R_plus = n_out_plus / (n_out_plus + n_out_minus)
    R_minus = n_out_minus / (n_out_plus + n_out_minus)

    # Calculate final operator probabilities
    if not test:
        P_out_plus = (P_SNR * R_plus) + ((1 - P_SNR) * 0.5)
        P_out_minus = (P_SNR * R_minus) + ((1 - P_SNR) * 0.5)
    else:
        P_out_plus_wo_multi = (P_SNR * R_plus) + ((1 - P_SNR) * 0.5)
        P_out_plus = 0.5 + (P_out_plus_wo_multi - 0.5) * np.exp(-(1 - eta) * n_tot)
        P_out_minus_wo_multi = (P_SNR * R_minus) + ((1 - P_SNR) * 0.5)
        P_out_minus = 0.5 + (P_out_minus_wo_multi - 0.5) * np.exp(-(1 - eta) * n_tot)

    return P_out_plus, P_out_minus
