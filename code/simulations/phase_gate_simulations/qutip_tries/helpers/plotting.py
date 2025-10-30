from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from typing import Callable, Dict, Any, Literal
from helpers.printing import print_data
from helpers.generic_computations import normalize_matrix
from helpers.compute_reflection_parameters import (
    compute_params,
    params_type,
    phase_shift,
)
from helpers.error_computation import (
    compute_reflection_amplitude_error,
)
from rich import print


def plot_qswitch_dynamics(
    tlist: np.ndarray[Any, np.dtype[np.float32]],
    result_0: qt.Result,
    result_1: qt.Result,
    a_out_0: np.ndarray,
    a_out_1: np.ndarray,
    real_input_shape: Callable[[np.ndarray, Dict[str, float]], np.ndarray],
    polarization: str,
    args: Dict[str, float],
) -> None:
    """
    Plot the dynamics of the quantum switch simulation.

    Parameters:
        tlist (np.ndarray): Time points.
        result_0 (Result): Simulation result for atom in |0⟩.
        result_1 (Result): Simulation result for atom in |1⟩.
        a_out_0 (np.ndarray): Output field for atom in |0⟩.
        a_out_1 (np.ndarray): Output field for atom in |1⟩.
        real_input_shape (Callable): Function for the input pulse envelope.
        args (dict): Arguments for the input pulse shape function.
    """
    fig, axs = plt.subplots(5, 2, figsize=(14, 16), sharex=True)

    # --- Population ---
    axs[0, 0].plot(tlist, result_0.expect[0], label="|0⟩")
    axs[0, 0].plot(tlist, result_0.expect[1], label="|1⟩")
    axs[0, 0].plot(tlist, result_0.expect[2], label="|e⟩")
    axs[0, 0].set_title("Atom in |0⟩: Populations")
    axs[0, 0].set_ylabel("Population")
    axs[0, 0].legend()
    axs[0, 0].grid()

    axs[0, 1].plot(tlist, result_1.expect[0], label="|0⟩")
    axs[0, 1].plot(tlist, result_1.expect[1], label="|1⟩")
    axs[0, 1].plot(tlist, result_1.expect[2], label="|e⟩")
    axs[0, 1].set_title("Atom in |1⟩: Populations")
    axs[0, 1].legend()
    axs[0, 1].grid()

    # --- Photon Number ---
    axs[1, 0].plot(tlist, result_0.expect[3], color="purple")
    axs[1, 0].set_title(
        "Atom in |0⟩: " + r"$\langle a^{\dagger}a\rangle = \langle n \rangle$"
    )
    axs[1, 0].set_ylabel("Photon Number")
    axs[1, 0].grid()

    axs[1, 1].plot(tlist, result_1.expect[3], color="purple")
    axs[1, 1].set_title(
        "Atom in |1⟩: " + r"$\langle a^{\dagger}a\rangle = \langle n \rangle$"
    )
    axs[1, 1].grid()

    # --- Cavity Field Re/Im ---
    axs[2, 0].plot(tlist, np.real(a_out_0), label="Re(a_out)", color="darkorange")
    axs[2, 0].plot(tlist, np.imag(a_out_0), label="Im(a_out)", color="teal")
    axs[2, 0].set_title("Atom in |0⟩: Output Field")
    axs[2, 0].set_ylabel("Field Amplitude")
    axs[2, 0].legend()
    axs[2, 0].grid()

    axs[2, 1].plot(tlist, np.real(a_out_1), label="Re(a_out)", color="darkorange")
    axs[2, 1].plot(tlist, np.imag(a_out_1), label="Im(a_out)", color="teal")
    axs[2, 1].set_title("Atom in |1⟩: Output Field")
    axs[2, 1].legend()
    axs[2, 1].grid()

    # --- Phase of Output ---
    axs[3, 0].plot(tlist, np.angle(a_out_0) / np.pi)
    axs[3, 0].set_title("Atom in |0⟩: Phase of Output Field")
    axs[3, 0].set_ylabel("Arg(a_out) / π")
    axs[3, 0].grid()

    axs[3, 1].plot(tlist, np.angle(a_out_1) / np.pi)
    axs[3, 1].set_title("Atom in |1⟩: Phase of Output Field")
    axs[3, 1].grid()

    # --- Input Envelope ---
    input_env = real_input_shape(tlist, args)
    axs[4, 0].plot(tlist, input_env, color="gray")
    axs[4, 0].set_title("Input Pulse Envelope")
    axs[4, 0].set_ylabel("Amplitude")
    axs[4, 0].set_xlabel("Time")
    axs[4, 0].grid()

    axs[4, 1].plot(tlist, input_env, color="gray")
    axs[4, 1].set_title("Input Pulse Envelope")
    axs[4, 1].set_xlabel("Time")
    axs[4, 1].grid()

    fig.suptitle(f"Polarization of Photon: {polarization}", fontsize=16)
    # plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


def plot_photon_number_and_population(
    tlist: np.ndarray, result_0: qt.Result, result_1: qt.Result
) -> None:
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], hspace=0.4, wspace=0.3)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.set_title(r"Driving only with $\pi$")
    ax1.plot(
        tlist, result_0.e_data["n_cav_pi"], linestyle="-", label=r"$n_{\pi,|0\rangle}$"
    )
    ax1.plot(
        tlist, result_1.e_data["n_cav_pi"], linestyle="--", label=r"$n_{\pi,|1\rangle}$"
    )
    ax1.plot(
        tlist, result_0.e_data["n_cav_v"], linestyle="-.", label=r"$n_{V,|0\rangle}$"
    )
    ax1.plot(
        tlist, result_1.e_data["n_cav_v"], linestyle=":", label=r"$n_{V,|1\rangle}$"
    )
    ax1.legend()
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Photon number")

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title(r"Population of $|0\rangle$")
    ax2.plot(tlist, result_0.e_data["P(0)"], linestyle="-", label=r"$|0\rangle$")
    ax2.plot(tlist, result_0.e_data["P(1)"], linestyle="--", label=r"$|1\rangle$")
    ax2.plot(tlist, result_0.e_data["P(e)"], linestyle="-.", label=r"$|e\rangle$")
    ax2.legend()
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Population")

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.set_title(r"Population of $|1\rangle$")
    ax3.plot(tlist, result_1.e_data["P(0)"], linestyle="--", label=r"$|0\rangle$")
    ax3.plot(tlist, result_1.e_data["P(1)"], linestyle="-", label=r"$|1\rangle$")
    ax3.plot(tlist, result_1.e_data["P(e)"], linestyle="-.", label=r"$|e\rangle$")
    ax3.legend()
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Population")

    plt.tight_layout()
    plt.show()


def styled_3d_bar(
    ax,
    mat,
    labels,
    cmap=LinearSegmentedColormap.from_list("my_list", ["white", "green"]),
    vmax=None,
    errors=None,
    cap_width=0.3,
):
    """
    Draw styled 3D bars on the given Axes3D `ax` based on `mat` (2D array).
    """
    mat = np.array(mat, dtype=float)
    n_rows, n_cols = mat.shape
    _x, _y = np.meshgrid(np.arange(n_cols), np.arange(n_rows))
    x = _x.ravel()
    y = _y.ravel()
    z = np.zeros_like(x)
    dz = mat.ravel()

    dx = dy = 0.8  # bar footprint

    # Choose colors: use a single-blue colormap but make zeros fully transparent
    if vmax is None:
        vmax = max(dz.max(), 1e-9)
    norm = plt.Normalize(0, vmax)
    cmap_obj = plt.get_cmap(cmap)
    colors = []
    for val in dz:
        if val <= 0:
            colors.append((1, 1, 1, 0.0))  # fully transparent
        else:
            rgba = cmap_obj(norm(val))
            # slightly darken edges: use same rgba but we will draw edges in black
            colors.append(rgba)

    # Draw bars: use edgecolor black for visibility and alpha from colormap
    ax.bar3d(
        x,
        y,
        z,
        dx,
        dy,
        dz,
        color=colors,
        edgecolor="gray",
        linewidth=0.6,
        shade=True,
        zsort="average",
        alpha=0.7,
    )
    if errors is not None:
        errors = np.array(errors, dtype=float).ravel()
        half_cap = cap_width / 2.0

        for xi, yi, zi, err in zip(x, y, dz, errors):
            if err > 0:
                top_z = zi
                err_low = max(top_z - err, 0)
                err_high = top_z + err

                # Vertical line
                ax.plot(
                    [xi + dx / 2.0] * 2,
                    [yi + dy / 2.0] * 2,
                    [err_low, err_high],
                    color="black",
                    linewidth=1.2,
                )

                # Optional cap (horizontal marker)
                ax.plot(
                    [xi + dx / 2.0 - half_cap, xi + dx / 2.0 + half_cap],
                    [yi + dy / 2.0, yi + dy / 2.0],
                    [err_high, err_high],
                    color="black",
                    linewidth=1.2,
                )

    # Style the base grid and panes to look clean
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    # ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("white")
    ax.yaxis.pane.set_edgecolor("white")
    # ax.zaxis.pane.set_edgecolor("white")

    # Ticks and labels
    ax.set_xticks(np.arange(n_cols) + dx / 2.0)
    ax.set_yticks(np.arange(n_rows) + dy / 2.0)
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Input", labelpad=15)
    ax.set_ylabel("Output", labelpad=15)
    ax.set_zlabel("Amplitude (a.u.)", labelpad=6)

    # Add numeric labels above nonzero bars (format to 2 decimals or 3 if small)
    if errors is not None:
        err_arr = np.array(errors, dtype=float).ravel()
    else:
        err_arr = np.zeros_like(dz)

    for xi, yi, zi, err in zip(x, y, dz, err_arr):
        if zi > 0:
            zpos = zi + vmax * 0.05  # slightly higher to avoid cap overlap

            # smart rounding
            fmt = "{:.2f}" if zi >= 0.01 else "{:.3f}"
            val_str = fmt.format(zi)
            err_str = fmt.format(err)

            if err > 0:
                label = f"{val_str} ± {err_str}"
            else:
                label = val_str

            ax.text(
                xi + dx / 2.0,
                yi + dy / 2.0,
                zpos,
                label,
                ha="center",
                va="bottom",
                fontsize=9,
                color="black",
                weight="bold",
            )  # Set limits and view
    ax.set_zlim(0, vmax * 1.15)
    # ax.view_init(elev=20, azim=45)


def plot_two_matrices_styled_grid(
    mat_left,
    mat_right,
    labels_left,
    labels_right,
    title_left,
    title_right,
    err_left: np.ndarray | None = None,
    err_right: np.ndarray | None = None,
    title_figure=None,
):
    """
    Plot two matrices (and their normalized versions) in a 2x2 grid.
    Left column = first matrix (raw + normalized)
    Right column = second matrix (raw + normalized)
    """
    # Auto-normalize if not provided
    norm_left = normalize_matrix(mat_left)
    norm_right = normalize_matrix(mat_right)

    # Create figure and axes (2x2)
    fig = plt.figure(figsize=(11, 8))
    axs = [
        fig.add_subplot(2, 2, 1, projection="3d"),
        fig.add_subplot(2, 2, 3, projection="3d"),
        fig.add_subplot(2, 2, 2, projection="3d"),
        fig.add_subplot(2, 2, 4, projection="3d"),
    ]

    # --- Left column: C-PHASE ---
    if err_left is not None:
        styled_3d_bar(axs[0], mat_left, labels_left, errors=err_left)
    else:
        styled_3d_bar(axs[0], mat_left, labels_left)
    axs[0].set_title(title_left, pad=10, fontsize=12)

    if err_right is not None:
        styled_3d_bar(axs[1], norm_left, labels_left, errors=err_right, vmax=1.0)
    else:
        styled_3d_bar(axs[1], norm_left, labels_left, vmax=1.0)
    axs[1].set_title(f"{title_left} - Normalized", pad=10, fontsize=12)

    # --- Right column: CNOT ---
    styled_3d_bar(axs[2], mat_right, labels_right)
    axs[2].set_title(title_right, pad=10, fontsize=12)

    styled_3d_bar(axs[3], norm_right, labels_right, vmax=1.0)
    axs[3].set_title(f"{title_right} - Normalized", pad=10, fontsize=12)

    # Adjust layout
    if title_figure is not None:
        plt.suptitle(title_figure, fontsize=16, y=0.98)
    # plt.tight_layout()
    plt.show()
    return norm_left, norm_right


def plot_cphase_and_cnot(
    basis: Dict[
        Literal["|0,pi>", "|1,pi>", "|0,V>", "|1,V>"],
        Dict[Literal["Delta a", "Delta c"], float],
    ],
    title: str,
    params_dir: params_type,
    params_dir_err: params_type,
    attenuate_light: bool = False,
):
    results = {}

    for label, detunings in basis.items():
        reflection_params = {
            "mu_rf": params_dir["mu_rf"],
            "mu_fc": params_dir["mu_fc"],
            "mu_fc_phi": params_dir["mu_fc_phi"],
            "kappa": params_dir["kappa"],
            "kappa_oc": params_dir["kappa_oc"],
            "d_w_r": detunings["Delta c"],
            "d_w_a": detunings["Delta a"],
            "gamma": params_dir["gamma"],
            "g": params_dir["g"],
        }
        reflection_amplitude, phase = compute_params(**reflection_params)
        reflection_amplitude_err = compute_reflection_amplitude_error(
            params_dir, params_dir_err, detunings["Delta a"], detunings["Delta c"]
        )
        reflection_phase_err = phase_shift(reflection_amplitude_err)
        power = np.abs(reflection_amplitude) ** 2
        results[label] = {
            "r": reflection_amplitude,
            "dr": reflection_amplitude_err,
            "|r|^2": power,
            "phase_rad": phase,
            "phase_deg": np.degrees(phase),
            "phase_rad_err": reflection_phase_err,
            "phase_deg_err": np.degrees(reflection_phase_err),
        }
    if attenuate_light:
        r_0_pi = results["|0,pi>"]["|r|^2"]
        r_1_pi = results["|1,pi>"]["|r|^2"]
        r_0_V = results["|0,V>"]["|r|^2"]
        r_1_V = results["|1,V>"]["|r|^2"]

        avg_pi = (r_0_pi + r_1_pi) / 2
        avg_v = (r_0_V + r_1_V) / 2
        reduction_v_polarization = avg_pi / avg_v

        results["|0,V>"]["r"] *= np.sqrt(reduction_v_polarization)
        results["|1,V>"]["r"] *= np.sqrt(reduction_v_polarization)
        results["|0,V>"]["dr"] *= np.sqrt(reduction_v_polarization)
        results["|1,V>"]["dr"] *= np.sqrt(reduction_v_polarization)
        results["|0,V>"]["|r|^2"] *= reduction_v_polarization
        results["|1,V>"]["|r|^2"] *= reduction_v_polarization

    c_0_pi = results["|0,pi>"]["r"]
    c_1_pi = results["|1,pi>"]["r"]
    c_0_V = results["|0,V>"]["r"]
    c_1_V = results["|1,V>"]["r"]

    d_0_pi = results["|0,pi>"]["dr"]
    d_1_pi = results["|1,pi>"]["dr"]
    d_0_V = results["|0,V>"]["dr"]
    d_1_V = results["|1,V>"]["dr"]

    r_0_pi = results["|0,pi>"]["|r|^2"]
    r_1_pi = results["|1,pi>"]["|r|^2"]
    r_0_V = results["|0,V>"]["|r|^2"]
    r_1_V = results["|1,V>"]["|r|^2"]

    dr_0_pi = np.abs(d_0_pi) ** 2
    dr_1_pi = np.abs(d_1_pi) ** 2
    dr_0_V = np.abs(d_0_V) ** 2
    dr_1_V = np.abs(d_1_V) ** 2

    print_data(title, basis, results)

    input_matrix = np.array(
        [
            [r_0_pi, 0.0, 0.0, 0.0],
            [0.0, r_1_pi, 0.0, 0.0],
            [0.0, 0.0, r_0_V, 0.0],
            [0.0, 0.0, 0.0, r_1_V],
        ]
    )
    input_matrix_err = np.array(
        [
            [dr_0_pi, 0.0, 0.0, 0.0],
            [0.0, dr_1_pi, 0.0, 0.0],
            [0.0, 0.0, dr_0_V, 0.0],
            [0.0, 0.0, 0.0, dr_1_V],
        ]
    )
    cnot_matrix = np.array(
        [
            [
                np.abs((c_1_pi + c_1_V) / 2) ** 2,
                np.abs((-c_1_pi + c_1_V) / 2) ** 2,
                0,
                0,
            ],
            [
                np.abs((-c_1_pi + c_1_V) / 2) ** 2,
                np.abs((c_1_pi + c_1_V) / 2) ** 2,
                0,
                0,
            ],
            [
                0,
                0,
                np.abs((+c_0_pi + c_0_V) / 2) ** 2,
                np.abs((-c_0_pi + c_0_V) / 2) ** 2,
            ],
            [
                0,
                0,
                np.abs((-c_0_pi + c_0_V) / 2) ** 2,
                np.abs((+c_0_pi + c_0_V) / 2) ** 2,
            ],
        ]
    )
    cnot_matrix_err = np.array(
        [
            [
                np.abs((d_1_pi + d_1_V) / 2) ** 2,
                np.abs((-d_1_pi + d_1_V) / 2) ** 2,
                0,
                0,
            ],
            [
                np.abs((-d_1_pi + d_1_V) / 2) ** 2,
                np.abs((d_1_pi + d_1_V) / 2) ** 2,
                0,
                0,
            ],
            [
                0,
                0,
                np.abs((+d_0_pi + d_0_V) / 2) ** 2,
                np.abs((-d_0_pi + d_0_V) / 2) ** 2,
            ],
            [
                0,
                0,
                np.abs((-d_0_pi + d_0_V) / 2) ** 2,
                np.abs((+d_0_pi + d_0_V) / 2) ** 2,
            ],
        ]
    )

    c_phase_labels = ["|0,π⟩", "|1,π⟩", "|0,V⟩", "|1,V⟩"]
    cnot_labels = ["|1,+⟩", "|1,-⟩", "|0,+⟩", "|0,-⟩"]

    normalized_cphase, normalized_cnot = plot_two_matrices_styled_grid(
        mat_left=input_matrix,
        err_left=input_matrix_err,
        labels_left=c_phase_labels,
        title_left="Reflection Intensity (C-PHASE)",
        mat_right=cnot_matrix,
        err_right=cnot_matrix_err,
        labels_right=cnot_labels,
        title_right="Population Probability (CNOT)",
        title_figure=title,
    )
    return normalized_cphase, normalized_cnot
