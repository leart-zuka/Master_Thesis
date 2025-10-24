import numpy as np


def normalize_matrix(input_matrix: np.ndarray):
    normalized = np.zeros_like(input_matrix)
    col_sums = input_matrix.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    normalized = input_matrix / col_sums
    return normalized
