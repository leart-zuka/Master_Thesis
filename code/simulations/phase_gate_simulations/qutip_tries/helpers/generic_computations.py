import numpy as np


def normalize_matrix(input_matrix: np.ndarray):
    """
    Note: Make sure that all entries are in terms of probabilities and not just raw complex numbers
          as that can lead to the normalization costant being close to 0, ergo blowing up your term
          which you're trying to normaliz :>
    """
    normalized = np.zeros_like(input_matrix)
    normalization_constants = normalizations(input_matrix)
    normalized = input_matrix / normalization_constants
    return normalized


def normalizations(input_matrix: np.ndarray):
    col_sums = input_matrix.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    return np.abs(col_sums)
