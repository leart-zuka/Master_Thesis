from typing import Dict, Tuple
import numpy as np
from helpers.printing import print_data
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D plotting)
from matplotlib import cm


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

input_matrix = np.array(
    [
        [c_0_pi, 0.0, 0.0, 0.0],
        [0.0, c_1_pi, 0.0, 0.0],
        [0.0, 0.0, c_0_V, 0.0],
        [0.0, 0.0, 0.0, c_1_V],
    ],
    dtype=np.float32,
)


def normalize_matrix(input_matrix: np.ndarray):
    normalized = np.zeros_like(input_matrix)
    for i, row in enumerate(input_matrix):
        normalization_constant = row.sum()
        normalized[i, :] = row / normalization_constant
    return normalized


normalized_input = normalize_matrix(input_matrix)


def plot_3d_power_matrix(matrix, labels=None, title_prefix="Reflection Power"):
    """
    Plots 3D bar plots for raw and row-normalized power matrices.

    Parameters:
        matrix : 2D np.array (square or rectangular)
            Matrix of power values (|r|^2)
        labels : list of str
            Labels for both rows and columns (length must match matrix.shape[0])
        title_prefix : str
            Prefix for plot titles
    """
    matrix = np.array(matrix, dtype=float)
    n_rows, n_cols = matrix.shape

    # Normalize rows safely
    normalized_matrix = np.zeros_like(matrix)
    for i, row in enumerate(matrix):
        s = row.sum()
        normalized_matrix[i, :] = row / s if s != 0 else 0

    # If no labels given, generate simple indices
    if labels is None:
        labels = [str(i) for i in range(n_rows)]

    # Flatten coordinates for bar3d
    _x, _y = np.meshgrid(np.arange(n_cols), np.arange(n_rows))
    x, y = _x.ravel(), _y.ravel()
    z = np.zeros_like(x)
    dx = dy = 0.6

    # Helper function to plot a single matrix
    def plot_single(mat, title):
        dz = mat.ravel()
        norm = plt.Normalize(dz.min() + 1e-9, dz.max())
        colors = cm.plasma(norm(dz))

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax.bar3d(x, y, z, dx, dy, dz, color=colors, shade=True, zsort="average")

        # Labels
        ax.set_xticks(np.arange(n_cols))
        ax.set_yticks(np.arange(n_rows))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=10)
        ax.set_yticklabels(labels, fontsize=10)
        ax.set_xlabel("Output Basis")
        ax.set_ylabel("Input Basis")
        ax.set_zlabel("Power (|r|²)")
        ax.set_title(title, fontsize=14, pad=12)

        # Annotate nonzero bars
        for xi, yi, zi in zip(x, y, dz):
            if zi > 0:
                ax.text(
                    xi,
                    yi,
                    zi + dz.max() * 0.03,
                    f"{zi:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    weight="bold",
                )

        ax.view_init(elev=25, azim=45)

        # Colorbar
        mappable = cm.ScalarMappable(norm=norm, cmap=cm.plasma)
        mappable.set_array(dz)
        fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.1, label="Power (|r|²)")

        plt.tight_layout()
        return fig, ax

    # Plot raw matrix
    plot_single(matrix, f"{title_prefix} — Raw")

    # Plot row-normalized matrix
    plot_single(normalized_matrix, f"{title_prefix} — Row-normalized")
    plt.show()


labels = ["|0,π⟩", "|1,π⟩", "|0,V⟩", "|1,V⟩"]
plot_3d_power_matrix(input_matrix, labels=labels, title_prefix="Reflection Power")
