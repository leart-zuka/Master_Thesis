from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rich import print


def overlap_poincare(
    folder,
    mode: Literal["_reflection", "_transmission", ""],
    version: str,
    s1_col=" Normalized s 1 ",
    s2_col=" Normalized s 2 ",
    s3_col=" Normalized s 3 ",
):
    states = ["H", "V", "A", "D", "R", "L"]
    # Define ideal target vectors for each input state
    ideals = {
        "H": [1, 0, 0],
        "V": [-1, 0, 0],
        "A": [0, -1, 0],
        "D": [0, 1, 0],
        "R": [0, 0, 1],
        "L": [0, 0, -1],
    }

    overlaps = []
    errors = []

    for state in states:
        # Read data
        path = f"./{folder}/{state}{mode}{version}.csv"
        df = pd.read_csv(path, sep=";", skiprows=8)

        # Calculate means (s) and standard errors (sigma)
        s_mean = df[[s1_col, s2_col, s3_col]].mean().values
        s_sem = df[[s1_col, s2_col, s3_col]].sem().values  # Standard Error of Mean

        # Get ideal vector (p)
        # Note: Your original code used s_d_ideal for 'A' and s_a_ideal for 'D', etc.
        # I have maintained that logic below.
        if state == "A":
            p = ideals["D"]
        elif state == "D":
            p = ideals["A"]
        elif state == "R":
            p = ideals["L"]
        elif state == "L":
            p = ideals["R"]
        else:
            p = ideals[state]

        # Calculate overlap and propagated error
        val, err = overlap_polarizations_with_error(s_mean, s_sem, p)

        overlaps.append(val)
        errors.append(err)
        print(f"Overlap {state}: {val:.4f} ± {err:.4f}")

    avg_overlap = np.mean(overlaps)
    # Error of the average (sum of squares / N)
    avg_error = np.sqrt(np.sum(np.square(errors))) / len(errors)

    print("-----------------------")
    print(f"Average Overlap: {avg_overlap:.4f} ± {avg_error:.4f}")
    print("-----------------------")


def overlap_polarizations(s, p):
    return 1 / 2 * (1 + np.dot(s, p))


def overlap_polarizations_with_error(s, sem, p):
    """
    Calculates overlap and the propagated standard error.
    s: mean vector [s1, s2, s3]
    sem: standard error vector [err_s1, err_s2, err_s3]
    p: ideal vector [p1, p2, p3]
    """
    overlap = 0.5 * (1 + np.dot(s, p))
    # Propagated error: 0.5 * sqrt( sum( (p_i * sigma_i)^2 ) )
    error = 0.5 * np.sqrt(np.sum((np.array(p) * np.array(sem)) ** 2))
    return overlap, error


print("6BW")
overlap_poincare(folder="6_BW/", mode="_reflection", version="")


print("7BW")
overlap_poincare(folder="7_BW/", mode="_reflection", version="")
