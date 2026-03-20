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
    r_dc=1e-4,
    eta=0.74205,
    pulse_length: float = 1e-6,
    detector_dead_time: float = 20e-9,
):

    n_out_plus = np.sum(np.abs(out_plus) ** 2) * dt
    n_out_minus = np.sum(np.abs(out_minus) ** 2) * dt

    l_plus = n_out_plus * eta
    l_minus = n_out_minus * eta

    gamma_p = l_plus + r_dc
    gamma_m = l_minus + r_dc

    # Probability that two events happen within the detector dead time
    prob_distinction = (
        2 * pulse_length * detector_dead_time - detector_dead_time**2
    ) / pulse_length**2

    # Probability that exactly 1 click was produced, consists of
    # two cases:
    # actually one event
    # + two events that happen within detector dead-time

    F_photons = (l_plus + l_minus) / (
        gamma_p * (1 + (gamma_p / 2) * prob_distinction)
        + gamma_m * (1 + (gamma_m / 2) * prob_distinction)
    )

    R_p = n_out_plus / (n_out_plus + n_out_minus)
    R_m = n_out_minus / (n_out_plus + n_out_minus)

    # Calculate final operator probabilities
    P_out_plus = F_photons * R_p + (1 - F_photons) * 0.5
    P_out_minus = F_photons * R_m + (1 - F_photons) * 0.5

    return P_out_plus, P_out_minus


def calculate_detection_probabilities_other(
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

    print(f"R_p: {R_p}")
    print(f"R_m: {R_m}")

    # Calculate final operator probabilities
    P_out_plus = F_photons * R_p + (1 - F_photons) * 0.5
    P_out_minus = F_photons * R_m + (1 - F_photons) * 0.5

    return P_out_plus, P_out_minus


# def calculate_detection_probabilities(
#     alpha_out_plus,  # coherent output field
#     alpha_out_minus,
#     n_quantum_plus,  # quantum-corrected photon flux
#     n_quantum_minus,
#     dt,
#     r_dc=1e-4,
#     eta=0.74205,
#     pulse_length=1e-6,
#     detector_dead_time=20e-9,
# ):
#     # Coherent (good) photon numbers
#     n_coherent_plus = np.sum(np.abs(alpha_out_plus) ** 2) * dt
#     n_coherent_minus = np.sum(np.abs(alpha_out_minus) ** 2) * dt
#
#     # Quantum excess (bad, from multi-photon saturation)
#     n_excess_plus = np.sum(n_quantum_plus) * dt - n_coherent_plus
#     n_excess_minus = np.sum(n_quantum_minus) * dt - n_coherent_minus
#
#     # Detected photon numbers
#     l_good_plus = n_coherent_plus * eta
#     l_good_minus = n_coherent_minus * eta
#     l_bad_plus = n_excess_plus * eta
#     l_bad_minus = n_excess_minus * eta
#
#     # Total events per detector (good + bad photons + dark counts)
#     gamma_p = (l_good_plus + l_bad_plus) + r_dc
#     gamma_m = (l_good_minus + l_bad_minus) + r_dc
#
#     prob_distinction = (
#         2 * pulse_length * detector_dead_time - detector_dead_time**2
#     ) / pulse_length**2
#
#     # F_good now only counts coherent (single-photon) events as "good"
#     # Both excess photons and dark counts contribute to "bad" events
#     F_photons = (l_good_plus + l_good_minus) / (
#         gamma_p * (1 + (gamma_p / 2) * prob_distinction)
#         + gamma_m * (1 + (gamma_m / 2) * prob_distinction)
#     )
#
#     # R± uses only the coherent part (the properly-gated photons)
#     n_coherent_tot = n_coherent_plus + n_coherent_minus
#     R_p = n_coherent_plus / n_coherent_tot if n_coherent_tot > 0 else 0.5
#     R_m = n_coherent_minus / n_coherent_tot if n_coherent_tot > 0 else 0.5
#
#     P_out_plus = F_photons * R_p + (1 - F_photons) * 0.5
#     P_out_minus = F_photons * R_m + (1 - F_photons) * 0.5
#     return P_out_plus, P_out_minus
