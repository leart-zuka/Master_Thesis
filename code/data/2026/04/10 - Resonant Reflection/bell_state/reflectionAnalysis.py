from __future__ import annotations

import os
import pickle
from typing import List
import numpy as np
import pandas as pd
from rich.progress import track
from rich import print
from rich.table import Table
from rich.console import Console
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit, minimize
from helper.fitting import R_coupled
from dataclasses import dataclass
from helper.numba_functions import (
    get_binary_up_and_down,
    get_data_from_main_h5_file,
    loop_over_time_stamps,
    get_atom_in_and_out_index,
    group_data_array,
)
from helper.plotting import channels_histo, plot_gate_3d
from helper.analysis_types import NormalModeSpectroscopyT, ReflectionGateT
from typing import Callable, Literal
from mpl_toolkits.mplot3d import Axes3D


_ATOM_KETS = {
    "Z": {
        "0": np.array([1, 0], dtype=complex),  # |0>
        "1": np.array([0, 1], dtype=complex),  # |1>
    },
    "X": {
        "0": np.array([1, 1], dtype=complex) / np.sqrt(2),  # |+>
        "1": np.array([1, -1], dtype=complex) / np.sqrt(2),  # |->
    },
    "Y": {
        "0": np.array([1, 1j], dtype=complex) / np.sqrt(2),  # |+i>
        "1": np.array([1, -1j], dtype=complex) / np.sqrt(2),  # |-i>
    },
}

# Photon kets in the {|V>, |pi>} basis (|V> is +1 of sigma_z^P by lab convention)
_PHOTON_KETS = {
    "HV": {
        "ch4": np.array([0, 1], dtype=complex),  # |pi>
        "ch7": np.array([1, 0], dtype=complex),  # |V>
    },
    "DA": {
        "ch4": np.array([1, 1], dtype=complex) / np.sqrt(2),  # |D>
        "ch7": np.array([1, -1], dtype=complex) / np.sqrt(2),  # |A>
    },
    "RL": {
        "ch4": np.array([1, 1j], dtype=complex) / np.sqrt(2),  # |R>
        "ch7": np.array([1, -1j], dtype=complex) / np.sqrt(2),  # |L>
    },
}


def _projector(setting: Setting, atom_outcome: str, channel: str) -> np.ndarray:
    """|psi_nu><psi_nu| for a given measurement outcome (atom, channel)."""
    atom_ket = _ATOM_KETS[setting.atom][atom_outcome]
    phot_ket = _PHOTON_KETS[setting.photon][channel]
    joint = np.kron(atom_ket, phot_ket)
    return np.outer(joint, joint.conj())


# ---------------------------------------------------------------------------
# T-matrix parametrization
# ---------------------------------------------------------------------------


def _t_from_params(t: np.ndarray) -> np.ndarray:
    """Build a 4x4 lower-triangular complex T from 16 real parameters."""
    T = np.zeros((4, 4), dtype=complex)
    T[0, 0] = t[0]
    T[1, 1] = t[1]
    T[2, 2] = t[2]
    T[3, 3] = t[3]
    T[1, 0] = t[4] + 1j * t[5]
    T[2, 0] = t[6] + 1j * t[7]
    T[2, 1] = t[8] + 1j * t[9]
    T[3, 0] = t[10] + 1j * t[11]
    T[3, 1] = t[12] + 1j * t[13]
    T[3, 2] = t[14] + 1j * t[15]
    return T


def _rho_from_params(t: np.ndarray) -> np.ndarray:
    """rho = T^dag T / Tr(T^dag T), automatically physical."""
    T = _t_from_params(t)
    G = T.conj().T @ T
    trace = np.trace(G).real
    if trace <= 0:
        return np.eye(4, dtype=complex) / 4
    return G / trace


# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------


def _collect_observations(
    probs_by_row: dict[int, Probs], inp: Input
) -> list[tuple[np.ndarray, float, float]]:
    """Gather (projector, p_measured, sigma) for each joint outcome."""
    obs: list[tuple[np.ndarray, float, float]] = []
    rows_for_input = [s for s in ALL_SETTINGS if s.inp == inp]
    for s in rows_for_input:
        pr = probs_by_row[s.row]
        for (atom_out, ch), (p, err) in pr.p.items():
            proj = _projector(s, atom_out, ch)
            sigma = max(err, 1e-6)
            obs.append((proj, p, sigma))
    return obs


def _neg_log_likelihood(
    t: np.ndarray, observations: list[tuple[np.ndarray, float, float]]
) -> float:
    """Gaussian NLL, dropping rho-independent constants."""
    rho = _rho_from_params(t)
    total = 0.0
    for proj, p_meas, sigma in observations:
        p_pred = np.real(np.trace(proj @ rho))
        total += (p_pred - p_meas) ** 2 / (2.0 * sigma**2)
    return total


# ---------------------------------------------------------------------------
# Top-level MLE routine
# ---------------------------------------------------------------------------


def mle_reconstruct(
    probs_by_row: dict[int, Probs],
    inp: Input,
    *,
    n_restarts: int = 3,
    rng: np.random.Generator | None = None,
    verbose: bool = False,
) -> tuple[np.ndarray, float]:
    """
    Find the physical rho maximizing the Gaussian likelihood of the data.

    Returns (rho, final_neg_log_likelihood).
    """
    if rng is None:
        rng = np.random.default_rng(0)

    observations = _collect_observations(probs_by_row, inp)

    # Initial guess: maximally mixed state (T = 0.5 * I)
    t0_identity = np.zeros(16)
    t0_identity[:4] = 0.5

    best_result = None
    for restart in range(n_restarts):
        t0 = (
            t0_identity if restart == 0 else t0_identity + 0.2 * rng.standard_normal(16)
        )
        result = minimize(
            _neg_log_likelihood,
            t0,
            args=(observations,),
            method="L-BFGS-B",
            options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-10},
        )
        if verbose:
            print(
                f"  restart {restart}: nll = {result.fun:.6f}, converged = {result.success}"
            )
        if best_result is None or result.fun < best_result.fun:
            best_result = result

    rho = _rho_from_params(best_result.x)
    rho = (rho + rho.conj().T) / 2  # numerical Hermiticity cleanup
    return rho, best_result.fun


# ---------------------------------------------------------------------------
# Driver with MC error bars
# ---------------------------------------------------------------------------


@dataclass
class MLEResult:
    rho: np.ndarray
    target: np.ndarray
    target_ket: np.ndarray
    fidelity_central: float
    fidelity_mc_mean: float
    fidelity_mc_std: float
    eigenvalues: np.ndarray
    purity: float
    neg_log_likelihood: float


def run_mle(
    probs_by_row: dict[int, Probs],
    *,
    mc_samples: int = 200,
    verbose: bool = True,
) -> dict[Input, MLEResult]:
    """Run MLE for both R and L with Monte Carlo error bars."""
    results: dict[Input, MLEResult] = {}
    rng = np.random.default_rng(42)

    for inp in ("R", "L"):
        if verbose:
            print(f"\n=== MLE reconstruction for {inp}-input ===")

        rho, nll = mle_reconstruct(probs_by_row, inp, verbose=verbose)

        psi = target_bell(inp)
        target_rho = np.outer(psi, psi.conj())
        F_central = fidelity(rho, psi)

        fids = []
        if mc_samples > 0:
            for sample_idx in range(mc_samples):
                sampled: dict[int, Probs] = {}
                for row, pr in probs_by_row.items():
                    raw = {
                        key: max(0.0, rng.normal(p, e)) for key, (p, e) in pr.p.items()
                    }
                    total = sum(raw.values())
                    if total <= 0:
                        sampled[row] = Probs(pr.p.copy())
                    else:
                        new_p = {
                            key: (raw[key] / total, e) for key, (p, e) in pr.p.items()
                        }
                        sampled[row] = Probs(new_p)
                try:
                    rho_s, _ = mle_reconstruct(
                        sampled, inp, n_restarts=1, verbose=False
                    )
                    fids.append(fidelity(rho_s, psi))
                except Exception:
                    continue
                if verbose and (sample_idx + 1) % 25 == 0:
                    print(
                        f"  MC sample {sample_idx + 1}/{mc_samples} (running F = {np.mean(fids):.4f})"
                    )

        F_mc_mean = float(np.mean(fids)) if fids else F_central
        F_mc_std = float(np.std(fids)) if fids else 0.0
        eigs = np.linalg.eigvalsh(rho).real
        purity = float(np.real(np.trace(rho @ rho)))

        results[inp] = MLEResult(
            rho=rho,
            target=target_rho,
            target_ket=psi,
            fidelity_central=F_central,
            fidelity_mc_mean=F_mc_mean,
            fidelity_mc_std=F_mc_std,
            eigenvalues=eigs,
            purity=purity,
            neg_log_likelihood=nll,
        )

        if verbose:
            print(f"\n  Density matrix (MLE):")
            print(np.array2string(rho, precision=3, suppress_small=True))
            print(f"  Eigenvalues: {np.round(eigs, 4)}")
            print(f"  Trace: {np.trace(rho).real:.4f}")
            print(f"  Purity: {purity:.4f}")
            print(f"  Fidelity (central): {F_central:.4f}")
            print(f"  Fidelity (MC): {F_mc_mean:.4f} +/- {F_mc_std:.4f}")
            print(f"  Final -log L: {nll:.4f}")

    return results


PhotonBasis = Literal["Vpi", "RL"]

_U_Vpi_to_RL = (1 / np.sqrt(2)) * np.array(
    [
        [1, 1],  # |V> component of |R>, |L>
        [1j, -1j],  # |pi> component of |R>, |L>
    ],
    dtype=complex,
)


def _change_photon_basis(rho_Vpi: np.ndarray) -> np.ndarray:
    """Transform 4x4 atom-photon rho from {V,pi} photon basis to {R,L}."""
    U_photon = _U_Vpi_to_RL  # 2x2
    U_full = np.kron(np.eye(2), U_photon)  # 4x4 acting on (atom)x(photon)
    return U_full.conj().T @ rho_Vpi @ U_full


def city_plot(
    rho: np.ndarray,
    *,
    photon_basis: PhotonBasis = "RL",
    target: np.ndarray | None = None,
    title: str = "",
    value: Literal["abs", "real", "imag"] = "abs",
    figsize: tuple[float, float] = (10, 10),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Render a 3D bar-chart city plot of a 4x4 density matrix.

    Parameters
    ----------
    rho : 4x4 complex density matrix in the {|0V>, |0pi>, |1V>, |1pi>} basis
          (as produced by density_matrix_from_stokes + project_to_physical).
    photon_basis : "Vpi" to plot as-is, "RL" to first rotate the photon to
          the R/L basis. Atom always stays in {|0>, |1>}.
    target : optional 4x4 ideal density matrix in the same basis; if given,
          drawn as open wireframe bars behind the solid data bars, Welte-style.
    value : which component of each matrix element to plot as bar height.
    """
    if photon_basis == "RL":
        rho_plot = _change_photon_basis(rho)
        labels = [r"$|0R\rangle$", r"$|0L\rangle$", r"$|1R\rangle$", r"$|1L\rangle$"]
        if target is not None:
            target_plot = _change_photon_basis(target)
    else:
        rho_plot = rho
        labels = [
            r"$|0V\rangle$",
            r"$|0\pi\rangle$",
            r"$|1V\rangle$",
            r"$|1\pi\rangle$",
        ]
        target_plot = target

    def _extract(m: np.ndarray) -> np.ndarray:
        if value == "abs":
            return np.abs(m)
        if value == "real":
            return m.real
        if value == "imag":
            return m.imag
        raise ValueError(f"Unknown value mode {value!r}")

    heights = _extract(rho_plot)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    # Bar positions
    _xpos, _ypos = np.meshgrid(np.arange(4), np.arange(4), indexing="ij")
    xpos = _xpos.flatten()
    ypos = _ypos.flatten()
    zpos = np.zeros_like(xpos, dtype=float)
    dx = dy = 0.7 * np.ones_like(xpos, dtype=float)
    dz = heights.flatten()

    # Solid data bars
    ax.bar3d(
        xpos, ypos, zpos, dx, dy, dz, color="crimson", alpha=0.8, edgecolor="black"
    )

    # Optional wireframe target bars
    if target_plot is not None:
        tgt_heights = _extract(target_plot).flatten()
        for xi, yi, hi in zip(xpos, ypos, tgt_heights):
            if abs(hi) < 1e-9:
                continue
            # Draw a rectangular wireframe at the target height
            _draw_wireframe_bar(ax, xi, yi, hi, dx[0], dy[0])

    ax.set_xticks(np.arange(4) + 0.35)
    ax.set_yticks(np.arange(4) + 0.35)
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)

    # ax.set_xlabel("Output", labelpad=20)
    # ax.set_ylabel("Input", labelpad=20)

    ax.set_zlim(0, max(0.50, dz.max() * 1.1))
    # ax.set_zlabel(
    #     {
    #         "abs": r"$|\rho_{ij}|$",
    #         "real": r"$\mathrm{Re}\,\rho_{ij}$",
    #         "imag": r"$\mathrm{Im}\,\rho_{ij}$",
    #     }[value]
    # )

    ax.set_zlabel("Absolute value", labelpad=15, rotation=90)

    # fig.tight_layout()
    # fig.subplots_adjust(left=0.05, right=0.85)
    fig.subplots_adjust(left=0.0, right=0.78, top=0.95, bottom=0.05)

    # if title:
    #     ax.set_title(title)
    ax.view_init(elev=15, azim=-30)
    if save_path:
        fig.savefig(save_path)

    return fig


def _draw_wireframe_bar(ax, x: float, y: float, h: float, dx: float, dy: float) -> None:
    """Draw the 12 edges of a rectangular bar as blue lines."""
    # Corners of the box footprint (counter-clockwise)
    corners_bottom = [
        (x, y, 0),
        (x + dx, y, 0),
        (x + dx, y + dy, 0),
        (x, y + dy, 0),
    ]
    corners_top = [(cx, cy, h) for (cx, cy, _) in corners_bottom]

    edges = []
    # Bottom
    for i in range(4):
        edges.append((corners_bottom[i], corners_bottom[(i + 1) % 4]))
    # Top
    for i in range(4):
        edges.append((corners_top[i], corners_top[(i + 1) % 4]))
    # Verticals
    for i in range(4):
        edges.append((corners_bottom[i], corners_top[i]))

    for p0, p1 in edges:
        xs = [p0[0], p1[0]]
        ys = [p0[1], p1[1]]
        zs = [p0[2], p1[2]]
        ax.plot(xs, ys, zs, color="blue", linewidth=1.0, alpha=0.6)


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------


def plot_bell_state_tomography(
    rho: np.ndarray,
    target: np.ndarray,
    inp: Literal["R", "L"],
    save_path: str | None = None,
) -> plt.Figure:
    """Thesis-ready default: |rho| in R/L photon basis with target wireframe."""
    return city_plot(
        rho,
        target=target,
        photon_basis="RL",
        value="abs",
        title=f"Atom-photon Bell state ({inp}-input)",
        save_path=save_path,
    )


class Analyzer:
    """
    Anal eyes
    Analyze
    Anal lies
    """

    def __init__(self, log_dir: str = "./", data_dir: str | None = "./") -> None:
        self.log_dir = log_dir
        self.data_dir = data_dir

        # ------ Definition of parameters ------
        self.sync_slow = 0
        self.sync_fast2 = 1  # QuTau Trigger (new trial within experiment run)
        self.lc_h = 2
        self.lc_v = 3
        self.kc_h = 4
        self.sync_fast = 5  # QuTau Trigger 3 (new experiment run)
        self.sd_trig = 6
        self.kc_v = 7

        self.ad_t = (
            0.13  # s - minimum atom trapping duration to be considered "good atom"
        )

        self.cooling_time = 25e-3
        self.fs_delay = 0.7e-6

        self.ps_save = True  # post selection save
        self.data_save = True

        # --- Colors ----
        self.colour = {
            "blueDark": (0, 0.3, 0.6),
            "blueLight": (0.5, 0.8, 1),
            "orangeDark": (1, 0.7, 0),
            "orangeLight": (1, 0.8, 0.6),
            "greenDark": (0, 0.6, 0.2),
            "greenLight": (0.7, 1, 0.5),
            "redDark": (0.9, 0, 0),
            "greyLight": (0.7, 0.7, 0.7),
        }

        self.data_dic = None
        self.data_arr = None

    def update_data_dir(self, new_data_dir: str) -> None:
        self.data_dir = new_data_dir

    def save_post_selection_data(self, base: Path, file_name: str, *args: pd.DataFrame):
        save_path = base / "goodAtomSelectorFiles"
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        with pd.ExcelWriter(f"{save_path}/{file_name}_atomParameters.xlsx") as writer:
            for arg in args:
                arg.to_excel(writer, sheet_name=f"{arg.sheet_name}")

    def post_selection(
        self,
        file_name: str,
        path: str | Path | None = None,
        file_type: str = ".h5",
        mean_kc_counts: int = 2000,
        no=10,
    ):
        if self.data_dir is None:
            raise Exception("Please define a data directory first")

        base = Path(path or self.data_dir)
        file_path = base / file_name

        # ------ We get the data ------
        self.data_arr = get_data_from_main_h5_file(base, file_name, file_type)
        wt_kc = 0.50 * mean_kc_counts  # witness threshold short cavity
        twot = 1.80 * mean_kc_counts  # two atom threshold
        atom_df = pd.DataFrame()
        data_photon_grouped = []
        data_time_grouped = []
        atom_in = []
        atom_out = []
        atom_in_histo = []
        atom_out_histo = []
        atoms_duration = []

        # ------ We enter the atom loop ------
        for atom_number in track(range(len(self.data_arr))):
            data_array = self.data_arr[atom_number]
            """
                    Really we just count all the all the counts in the short cavity for a run, and then we save:
                        counts in KC -> dataPhotonKC
                        times for an atom in KC -> dataTimeKC

                    timeStamps are all the time stamps in seconds
            """
            # --- Short Cavity --- #
            time_stamps: np.ndarray = data_array[self.sync_fast][1:-1]
            data_photon_kc, data_time_kc = loop_over_time_stamps(
                data_array, time_stamps
            )

            current_data_photon_grouped = group_data_array(data_photon_kc, no)
            current_data_time_grouped = group_data_array(data_time_kc, no)

            data_photon_grouped = np.concatenate(
                [data_photon_grouped, current_data_photon_grouped]
            )
            data_time_grouped = np.concatenate(
                [data_time_grouped, current_data_time_grouped]
            )

            atom_in_index, atom_out_index = get_atom_in_and_out_index(
                current_data_photon_grouped, wt_kc, twot
            )

            try:
                in_val = (
                    current_data_time_grouped[atom_in_index]
                    - data_array[self.sync_slow][0]
                )
                in_histo_val = current_data_time_grouped[atom_in_index]
            except Exception:
                in_val = 0
                in_histo_val = data_array[self.sync_slow][0]

            atom_in.append(in_val)
            atom_in_histo.append(in_histo_val)

            try:
                out_val = (
                    current_data_time_grouped[atom_out_index]
                    - data_array[self.sync_slow][0]
                )
                out_histo_val = current_data_time_grouped[atom_out_index]
            except Exception:
                out_val = 0
                out_histo_val = data_array[self.sync_slow][0]

            atom_out.append(out_val)
            atom_out_histo.append(out_histo_val)

            atoms_duration.append(atom_out[-1] - atom_in[-1])

        # %% - DATA ALLOCATION IN A DATA FRAME
        atom_df["atomsDuration"] = atoms_duration
        atom_df["atomsIn"] = atom_in
        atom_df["atomsOut"] = atom_out
        atom_df.sheet_name = "atomParameters"

        # Good atoms are selected, added in the data frame and in a dictionary
        """
            Only select the ones where the duration inside of the cavity is
            above a certain threshold
        """
        good_atoms_df = atom_df[(atom_df["atomsDuration"] >= self.ad_t)]
        good_atoms_df.sheet_name = "goodAtoms"
        good_atoms_dict = {
            i: [good_atoms_df["atomsIn"][i], good_atoms_df["atomsOut"][i]]
            for i in list(good_atoms_df.index)
        }
        good_atoms_dict_df = pd.DataFrame.from_dict(good_atoms_dict)
        good_atoms_dict_df.sheet_name = "goodAtomsDic"

        # The conditions for good atoms selection are saved in a data frame
        conds_df = pd.DataFrame()
        conds_df["Conditions"] = ["Single atom time threshold (s)"]
        conds_df["Bounds"] = [self.ad_t]
        conds_df.sheet_name = "gootAtomsConds"

        # %% ------ We plot the data ------
        plt.close("all")

        f = plt.figure(file_name, figsize=[17, 14])

        # --- kc counts plot --- #
        plt.plot(
            data_time_grouped,
            data_photon_grouped,
            color="tab:orange",
            label="Short Cavity counts",
            ls="None",
            marker=".",
        )
        plt.vlines(
            atom_in_histo, -20, 0, color="grey", linestyle="--", label="atom start time"
        )
        plt.vlines(
            atom_out_histo, -20, 0, color="red", linestyle="--", label="atom out time"
        )
        plt.hlines(
            wt_kc, atom_in_histo[0], atom_out_histo[-1], color="tab:green", alpha=0.2
        )
        plt.hlines(
            twot, atom_in_histo[0], atom_out_histo[-1], color="tab:red", alpha=0.2
        )

        for atom_number in range(len(atom_in_histo)):
            """
                if an atom is in the cavity and lives long enough the
                background will be dyed in a color in specific color
                """
            if atom_out_histo[atom_number] - atom_in_histo[atom_number] >= self.ad_t:
                plt.axvspan(
                    atom_in_histo[atom_number],
                    atom_out_histo[atom_number],
                    alpha=0.5,
                    color="tab:purple",
                )

            # print number of atom below
            plt.text(
                atom_in_histo[atom_number],
                -20 + 20 * (atom_number % 2),
                str(atom_number),
                fontsize=10,
            )

        plt.xlim(xmin=atom_in_histo[0] - 2, xmax=atom_out_histo[-1] + 2)
        plt.ylim(-1 * wt_kc, twot * 2)
        plt.legend()

        plt.tight_layout()
        # plt.show()

        # %% - DATA SAVING
        if self.ps_save is True:
            self.save_post_selection_data(
                base, file_name, atom_df, good_atoms_dict_df, good_atoms_df, conds_df
            )
            with open(
                f"{base}/goodAtomSelectorFiles/{file_name}_goodAtoms.pkl", "wb"
            ) as file:
                pickle.dump(good_atoms_dict, file)
            f.savefig(f"{base}/goodAtomSelectorFiles/{file_path}.png")

        return good_atoms_dict, atom_in_histo, atom_out_histo

    def get_trap_times(self, goodAtomsDic, atomInHisto, atomOutHisto):
        list_trappingDuration = []
        for key in goodAtomsDic:
            list_trappingDuration.append(goodAtomsDic[key][1] - goodAtomsDic[key][0])

        averageTrapTime = np.mean(list_trappingDuration)
        averageTrapTime_err = np.std(list_trappingDuration) / np.sqrt(
            np.size(list_trappingDuration)
        )
        trappingProbability = len(list_trappingDuration) / len(atomInHisto) * 100
        dutyCycle = (
            sum(list_trappingDuration) / (atomOutHisto[-1] - atomInHisto[0]) * 100
        )

        table = Table(title="Trapping Info")
        table.add_column("Attribute", justify="left", style="yellow")
        table.add_column("Value", justify="right", style="blue")

        table.add_row(
            "Average single atom trapping time",
            f"{averageTrapTime:.2f} +/- {averageTrapTime_err:.2f}",
        )
        table.add_row("Atom trapping probability", f"{trappingProbability:.0f}%")
        table.add_row("Duty cycle", f"{dutyCycle:.0f}%")
        console = Console()
        console.print(table)
        return averageTrapTime, averageTrapTime_err, trappingProbability, dutyCycle

    def normal_mode_spectroscopy(
        self,
        file_name: str,
        parameters: NormalModeSpectroscopyT,
        path: str | Path | None = None,
        file_type: str = ".h5",
        plot_histogramm: bool = False,
        fit_function: bool = False,
    ):
        if self.data_dir is None:
            raise Exception("please define a data directory fist")
        base = Path(path or self.data_dir)
        file_post_selected = (
            base / "goodAtomSelectorFiles" / f"{file_name}_goodAtoms.pkl"
        )
        if not os.path.exists(file_post_selected):
            self.post_selection(file_name, base, file_type)

        print(
            f"Analyzing Normal Mode Spectroscopy of file [green]{file_post_selected}[/green]"
        )

        with open(file_post_selected, "rb") as file:
            atom_dict = pickle.load(file)

        data_arr = self.data_good_atoms(atom_dict, base, file_name, file_type)

        # --- Load parameters ---
        trigger_delay = parameters["trigger_delay"]
        cooling_duration = parameters["cooling_duration"]
        optical_pumping_duration = parameters["optical_pumping_duration"]
        pulse_delay = parameters["pulse_delay_SD"]
        pulse_duration = parameters["pulse_duration_SD"]
        sequence_duration = parameters["sequence_duration"]
        frequency_span = parameters["frequency_span"]
        points_per_scan = parameters["points_per_scan"]
        trials_per_point = parameters["trials_per_point"]
        frequency_center = parameters["frequency_center"]

        binsize = 20 * 1e-9

        # sequence gates
        optical_pumping_gate = [
            trigger_delay + cooling_duration,
            trigger_delay + cooling_duration + optical_pumping_duration,
        ]
        write_gate = [
            optical_pumping_gate[1] + pulse_delay,
            optical_pumping_gate[1] + pulse_delay + pulse_duration,
        ]

        # plotting
        if plot_histogramm:
            binNum = int(sequence_duration / binsize)
            detectors = [self.kc_h, self.kc_v, self.lc_h, self.lc_v, self.sd_trig]
            colors = ["violet", "violet", "tab:blue", "tab:blue", "orange"]
            fsdelay = {
                self.kc_h: 0,
                self.kc_v: 12e-9,
                self.lc_h: 0,
                self.lc_v: 0.0,
                self.sd_trig: 0.0,
            }
            channels_histo(
                data_arr,
                detectors,
                write_gate,
                binNum,
                self.sync_fast,
                sequence_duration,
                fsdelay,
                file_name,
                colors,
            )
            plt.show(block=True)

        fast_sequence_step_size = (
            data_arr[self.sync_fast][-1] - data_arr[self.sync_fast][-2]
        )
        scan_duration = fast_sequence_step_size * points_per_scan * trials_per_point

        print(":waffle: Looping over [purple]Fast Sequence Triggers[/purple] now")

        binary_up, binary_down = get_binary_up_and_down(
            points_per_scan,
            trials_per_point,
            scan_duration,
            data_arr[self.kc_h],
            data_arr[self.kc_v],
            data_arr[self.sync_fast],
            data_arr[self.sync_fast2],
            write_gate,
        )

        binary_up_mean = [np.mean(point) for point in binary_up]
        binary_up_err = [np.std(point) / np.sqrt(len(point)) for point in binary_up]
        binary_down_mean = [np.mean(point) for point in binary_down]
        binary_down_err = [np.std(point) / np.sqrt(len(point)) for point in binary_down]
        frequency_span = (
            np.linspace(0, frequency_span, int(points_per_scan / 2)) - frequency_center
        )

        if fit_function:
            p0 = [
                # A (normalization, negative guess here to match scaling of data)
                -27,
                0.187,  # f_res (MHz offset of resonance)
                30,  # g (coupling strength, MHz)
                58,  # kappa (total cavity decay rate, MHz)
                58 * 0.85,  # kappa_oc (outcoupling, ~85% of total kappa)
                0.978,  # MM_rf (close to ideal)
                0.873,  # MM_fc (80% coupling in)
                3.0333,  # gamma (free space decay rate, MHz)
                0.01,  # offset (background level)
                0.0,  # a (slope term for detuning-dependent broadening)
            ]

            bounds = (
                [
                    -np.inf,
                    0.01,
                    10,
                    58,
                    49,
                    0.972,
                    0.871,
                    3.0318,
                    -np.inf,
                    -np.inf,
                ],
                [
                    np.inf,
                    0.3,
                    50,
                    59,
                    50,
                    0.984,
                    0.875,
                    3.0354,
                    np.inf,
                    np.inf,
                ],
            )
            popt, pcov = curve_fit(
                R_coupled,
                frequency_span,
                binary_up_mean + np.flip(binary_down_mean),
                p0,
                bounds=bounds,
                maxfev=10000,
            )

            pcov = np.sqrt(np.diag(pcov))
            plt.plot(
                frequency_span,
                R_coupled(frequency_span, *popt),
                label="Model fit",
                color="red",
                linewidth=3,
                linestyle="-.",
            )
            plt.suptitle(
                "\n $\\mathbf{g:\\ %.1f\\ MHz\\ \\pm\\ %.1f\\ MHz}$"
                % (
                    popt[2],
                    pcov[2],
                )
            )

        plt.errorbar(
            frequency_span,
            binary_up_mean + np.flip(binary_down_mean),
            binary_up_err + np.flip(binary_down_err),
            linestyle="",
            marker="o",
            label="Measurement data",
        )
        plt.show()

    def reflection_analysis(
        self,
        file_name: str,
        parameters: ReflectionGateT,
        path: str | Path | None = None,
        file_type: str = ".h5",
        post_select_sd: bool = False,
        plot_histogram: bool = False,
    ):
        if self.data_dir is None:
            raise Exception("please define a data directory first")

        base = Path(path or self.data_dir)

        file_post_selected = (
            base / "goodAtomSelectorFiles" / f"{file_name}_goodAtoms.pkl"
        )

        if not os.path.exists(file_post_selected):
            self.post_selection(file_name, base, file_type)

        print(f"Analyzing CNOT Gate of file [green]{file_post_selected}[/green]")

        with open(file_post_selected, "rb") as file:
            atom_dict = pickle.load(file)
        data_arr = self.data_good_atoms(atom_dict, base, file_name, file_type)

        trigger_delay = parameters["trigger_delay"]
        cooling_duration = parameters["cooling_duration"]
        optical_pumping_duration = parameters["optical_pumping_duration"]
        pulse_delay = parameters["pulse_delay"]
        pulse_duration = parameters["pulse_duration"]
        pulse_delay_SD_after = parameters["pulse_delay_SD_after"]
        pulse_duration_SD_after = parameters["pulse_duration_SD_after"]
        # pulse_delay_SD_before = parameters["pulse_delay_SD_before"]
        # pulse_duration_SD_before = parameters["pulse_duration_SD_before"]
        sequence_duration = parameters["sequence_duration"]

        binsize = 20 * 1e-9

        # sequence gates
        optical_pumping_gate = [
            trigger_delay + cooling_duration,
            trigger_delay + cooling_duration + optical_pumping_duration,
        ]
        write_gate = [
            optical_pumping_gate[1] + pulse_delay,  # Start of photon pulse
            optical_pumping_gate[1]
            + pulse_delay
            + pulse_duration,  # End of photon pulse
            optical_pumping_gate[1]
            + pulse_delay
            + pulse_duration
            + pulse_delay_SD_after,  # Start of SD
            optical_pumping_gate[1]
            + pulse_delay
            + pulse_duration
            + pulse_delay_SD_after
            + pulse_duration_SD_after,  # End of SD
            # optical_pumping_gate[1] + pulse_delay_SD_before,  # Start of photon pulse
            # optical_pumping_gate[1]
            # + pulse_delay_SD_before
            # + pulse_duration_SD_before
        ]

        # plotting
        if plot_histogram:
            binNum = int(sequence_duration / binsize)
            detectors = [self.kc_h, self.kc_v, self.lc_h, self.lc_v, self.sd_trig]
            colors = ["violet", "violet", "tab:blue", "tab:blue", "orange"]
            fsdelay = {
                self.kc_h: 0,
                self.kc_v: 12e-9,
                self.lc_h: 0,
                self.lc_v: 0.0,
                self.sd_trig: 0.0,
            }
            channels_histo(
                data_arr,
                detectors,
                write_gate,
                binNum,
                self.sync_fast,
                sequence_duration,
                fsdelay,
                file_name,
                colors,
            )
            plt.show(block=True)

        counts_ch_4_atom_1 = []
        counts_ch_7_atom_1 = []
        counts_ch_4_atom_0 = []
        counts_ch_7_atom_0 = []
        click_trials = 0

        # sync fast   is qutau trigger   (aka new trial)
        # sync fast 2 is qutau trigger 3 (aka new experiment)

        length = int(np.floor(len(data_arr[self.sync_fast][:-1]) / 2))
        # length = int(np.floor(len(data_arr[self.sync_fast][:-1])))
        print(length)

        for i in track(range(length)):
            t0 = data_arr[self.sync_fast][i]

            start = t0 + write_gate[0]
            end = t0 + write_gate[1]
            start_sd_after = t0 + write_gate[2]
            end_sd_after = t0 + write_gate[3]
            # start_sd_before = t0 + write_gate[4]
            # end_sd_before = t0 + write_gate[5]

            # start_4, end_4, start_4_sd_after, end_4_sd_after,start_4_sd_before, end_4_sd_before = np.searchsorted(
            #     data_arr[self.kc_h],
            #     [start, end, start_sd_after, end_sd_after, start_sd_before, end_sd_before]
            # )

            start_4, end_4, start_4_sd_after, end_4_sd_after = np.searchsorted(
                data_arr[self.kc_h], [start, end, start_sd_after, end_sd_after]
            )

            start_7, end_7, start_7_sd_after, end_7_sd_after = np.searchsorted(
                data_arr[self.kc_v],
                [start, end, start_sd_after, end_sd_after],
            )

            n4 = end_4 - start_4
            n7 = end_7 - start_7

            n4_sd = end_4_sd_after - start_4_sd_after
            n7_sd = end_7_sd_after - start_7_sd_after

            # ps4 = end_4_sd_before - start_4_sd_before
            # ps7 = end_7_sd_before - start_7_sd_before

            n_sd = n4_sd + n7_sd
            # ps = ps4 + ps7

            if n4 + n7 == 0 or n4 + n7 > 1:
                continue

            click_trials += 1

            if n_sd > 0:
                counts_ch_4_atom_1.append(n4)
                counts_ch_7_atom_1.append(n7)
            else:
                # pass
                counts_ch_4_atom_0.append(n4)
                counts_ch_7_atom_0.append(n7)

        print(click_trials)
        tot_clicks_4_atom_1 = np.sum(counts_ch_4_atom_1)
        tot_clicks_7_atom_1 = np.sum(counts_ch_7_atom_1)
        tot_clicks_4_atom_0 = np.sum(counts_ch_4_atom_0)
        tot_clicks_7_atom_0 = np.sum(counts_ch_7_atom_0)
        tot_clicks = (
            tot_clicks_4_atom_1
            + tot_clicks_7_atom_1
            + tot_clicks_4_atom_0
            + tot_clicks_7_atom_0
        )

        P_4_atom_0 = tot_clicks_4_atom_0 / tot_clicks
        P_4_atom_1 = tot_clicks_4_atom_1 / tot_clicks
        P_7_atom_0 = tot_clicks_7_atom_0 / tot_clicks
        P_7_atom_1 = tot_clicks_7_atom_1 / tot_clicks

        P_4_atom_0_err = np.sqrt(P_4_atom_0 * (1 - P_4_atom_0) / click_trials)
        P_4_atom_1_err = np.sqrt(P_4_atom_1 * (1 - P_4_atom_1) / click_trials)
        P_7_atom_0_err = np.sqrt(P_7_atom_0 * (1 - P_7_atom_0) / click_trials)
        P_7_atom_1_err = np.sqrt(P_7_atom_1 * (1 - P_7_atom_1) / click_trials)

        efficiency = click_trials / length
        try:
            efficiency_err = np.sqrt(efficiency * (1 - efficiency) / click_trials)
        except Exception as e:
            efficiency_err = 0

        self.data_arr = None

        return (
            P_4_atom_0,
            P_4_atom_0_err,
            P_4_atom_1,
            P_4_atom_1_err,
            P_7_atom_0,
            P_7_atom_0_err,
            P_7_atom_1,
            P_7_atom_1_err,
            efficiency,
            efficiency_err,
        )

    def data_good_atoms(
        self, atom_dict, base, file_name, file_type
    ) -> List[np.ndarray]:
        data_vd = [[] for _ in range(8)]
        if self.data_arr is None:
            self.data_arr = get_data_from_main_h5_file(base, file_name, file_type)

        for atom_number in atom_dict.keys():
            data_array = self.data_arr[atom_number]

            start_sync_fast_trigger = np.searchsorted(
                data_array[self.sync_fast],
                data_array[self.sync_slow][0] + atom_dict[atom_number][0],
            )

            end_sync_fast_trigger = np.searchsorted(
                data_array[self.sync_fast],
                data_array[self.sync_slow][0] + atom_dict[atom_number][1],
            )

            time_start = data_array[self.sync_fast][start_sync_fast_trigger] - 1e-9
            time_end = data_array[self.sync_fast][end_sync_fast_trigger] - 1e-9

            for i in range(8):
                start = np.searchsorted(data_array[i], time_start)
                end = np.searchsorted(data_array[i], time_end)

                data_vd[i] = np.append(data_vd[i], data_array[i][start:end])

        return data_vd


# ---------------------------------------------------------------------------
# Pauli matrices (standard conventions; the photon-Z sign flip is handled
# in the eigenvalue lookup, not in the matrices themselves)
# ---------------------------------------------------------------------------

_I2 = np.array([[1, 0], [0, 1]], dtype=complex)
_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULI = (_I2, _SX, _SY, _SZ)


# ---------------------------------------------------------------------------
# Channel <-> eigenvalue mapping
# ---------------------------------------------------------------------------

# eigenvalue of the photon Pauli observable for each (basis, channel).
# This is the ONE place to edit if your detector/waveplate conventions
# differ or a basis is physically flipped in hardware.
PhotonBasis = Literal["RL", "HV", "DA"]

PHOTON_EIGENVALUE: dict[PhotonBasis, dict[str, int]] = {
    # sigma_y: |R> -> +1, |L> -> -1
    "RL": {"ch4": +1, "ch7": -1},
    # sigma_z (lab-flipped): |V> -> +1, |pi> = |H> -> -1
    # If ch4 is the pi/H detector: ch4 -> -1. Adjust here if opposite.
    "HV": {"ch4": -1, "ch7": +1},
    # sigma_x: |D> -> +1, |A> -> -1
    "DA": {"ch4": +1, "ch7": -1},
}

# atom Pauli eigenvalue is the same in every basis (X, Y, Z are all measured
# via SD after an appropriate MW rotation, so the "raw" SD label keeps the
# same +/-1 assignment):
#   |0> -> +1,  |1> -> -1
ATOM_EIGENVALUE = {"0": +1, "1": -1}


# ---------------------------------------------------------------------------
# Measurement row bookkeeping
# ---------------------------------------------------------------------------

AtomBasis = Literal["Z", "X", "Y"]
Input = Literal["R", "L"]


@dataclass(frozen=True)
class Setting:
    """One of the 18 tomography settings."""

    row: int
    atom: AtomBasis
    photon: PhotonBasis
    inp: Input


def _build_row_map() -> list[Setting]:
    """Construct the full 18-row schedule matching the Hermes FPGA layout."""
    photon_blocks: list[PhotonBasis] = ["RL", "HV", "DA"]
    atom_cycle: list[AtomBasis] = ["Z", "X", "Y"]  # FPGA 0, 1, 2

    settings: list[Setting] = []
    row = 1
    for photon in photon_blocks:
        for atom in atom_cycle:
            for inp in ("R", "L"):
                settings.append(Setting(row=row, atom=atom, photon=photon, inp=inp))
                row += 1
    return settings


ALL_SETTINGS: list[Setting] = _build_row_map()


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------


@dataclass
class Probs:
    """Joint probabilities for a single measurement setting.

    Keys are (atom_outcome, channel). Values are (probability, uncertainty).
    All probabilities are in [0, 1] (not percentages).
    """

    p: dict[tuple[str, str], tuple[float, float]]

    @classmethod
    def from_reflection(
        cls,
        P_4_atom_0: float,
        P_4_atom_0_err: float,
        P_4_atom_1: float,
        P_4_atom_1_err: float,
        P_7_atom_0: float,
        P_7_atom_0_err: float,
        P_7_atom_1: float,
        P_7_atom_1_err: float,
    ) -> "Probs":
        """Wrap the 8 return values of analyzer.reflection_analysis()."""
        return cls(
            {
                ("0", "ch4"): (P_4_atom_0, P_4_atom_0_err),
                ("1", "ch4"): (P_4_atom_1, P_4_atom_1_err),
                ("0", "ch7"): (P_7_atom_0, P_7_atom_0_err),
                ("1", "ch7"): (P_7_atom_1, P_7_atom_1_err),
            }
        )

    def atom_marginal(self) -> dict[str, tuple[float, float]]:
        """P(atom=0), P(atom=1), summed over photon channels."""
        out: dict[str, tuple[float, float]] = {}
        for atom in ("0", "1"):
            p = sum(self.p[(atom, ch)][0] for ch in ("ch4", "ch7"))
            e = np.sqrt(sum(self.p[(atom, ch)][1] ** 2 for ch in ("ch4", "ch7")))
            out[atom] = (p, e)
        return out

    def photon_marginal(self) -> dict[str, tuple[float, float]]:
        """P(ch4), P(ch7), summed over atom outcomes."""
        out: dict[str, tuple[float, float]] = {}
        for ch in ("ch4", "ch7"):
            p = sum(self.p[(atom, ch)][0] for atom in ("0", "1"))
            e = np.sqrt(sum(self.p[(atom, ch)][1] ** 2 for atom in ("0", "1")))
            out[ch] = (p, e)
        return out


# ---------------------------------------------------------------------------
# Stokes parameter computation
# ---------------------------------------------------------------------------

# Pauli index convention: 0 = I, 1 = X, 2 = Y, 3 = Z
_ATOM_BASIS_TO_INDEX = {"X": 1, "Y": 2, "Z": 3}
_PHOTON_BASIS_TO_INDEX = {"DA": 1, "RL": 2, "HV": 3}


def _joint_correlator(probs: Probs, photon_basis: PhotonBasis) -> tuple[float, float]:
    """
    <sigma_i^A (x) sigma_j^P> for atom non-identity.

    Sum over 4 joint outcomes weighted by (atom_eig)(photon_eig).
    """
    eig_phot = PHOTON_EIGENVALUE[photon_basis]
    val = 0.0
    var = 0.0
    for atom in ("0", "1"):
        a_sign = ATOM_EIGENVALUE[atom]
        for ch in ("ch4", "ch7"):
            p_sign = eig_phot[ch]
            p, err = probs.p[(atom, ch)]
            val += a_sign * p_sign * p
            var += err**2  # coefficients are +/-1, square to +1
    return val, float(np.sqrt(var))


def _atom_marginal_correlator(probs: Probs) -> tuple[float, float]:
    """<sigma_i^A (x) I^P> = P(atom=0) - P(atom=1)."""
    m = probs.atom_marginal()
    val = m["0"][0] - m["1"][0]
    err = float(np.sqrt(m["0"][1] ** 2 + m["1"][1] ** 2))
    return val, err


def _photon_marginal_correlator(
    probs: Probs, photon_basis: PhotonBasis
) -> tuple[float, float]:
    """<I^A (x) sigma_j^P> = +1 channel prob - -1 channel prob."""
    eig = PHOTON_EIGENVALUE[photon_basis]
    m = probs.photon_marginal()
    val = sum(eig[ch] * m[ch][0] for ch in ("ch4", "ch7"))
    err = float(np.sqrt(sum(m[ch][1] ** 2 for ch in ("ch4", "ch7"))))
    return val, err


def _pool(values: list[tuple[float, float]]) -> tuple[float, float]:
    """Inverse-variance weighted mean of independent estimates."""
    weights = np.array([1.0 / e**2 if e > 0 else 0.0 for _, e in values])
    vals = np.array([v for v, _ in values])
    if weights.sum() == 0:
        return float(np.mean(vals)), float("inf")
    mean = float(np.sum(weights * vals) / weights.sum())
    err = float(np.sqrt(1.0 / weights.sum()))
    return mean, err


def compute_stokes(
    probs_by_row: dict[int, Probs],
    inp: Input,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Assemble the 4x4 S-matrix (and its errors) for one input polarization.

    S[i, j] corresponds to <sigma_i^A (x) sigma_j^P> with
      i,j = 0 -> I, 1 -> X, 2 -> Y, 3 -> Z.
    """
    S = np.zeros((4, 4))
    S_err = np.zeros((4, 4))

    # S00 = 1 exactly
    S[0, 0] = 1.0

    # Select rows for this input polarization
    rows_for_input = [s for s in ALL_SETTINGS if s.inp == inp]

    # --- Joint correlators S_ij (i, j in {1,2,3}) ---
    for s in rows_for_input:
        i = _ATOM_BASIS_TO_INDEX[s.atom]
        j = _PHOTON_BASIS_TO_INDEX[s.photon]
        val, err = _joint_correlator(probs_by_row[s.row], s.photon)
        S[i, j] = val
        S_err[i, j] = err

    # --- Atom marginals S_i0 (i in {1,2,3}): pool rows with matching atom basis ---
    for atom, i in _ATOM_BASIS_TO_INDEX.items():
        candidates = [s for s in rows_for_input if s.atom == atom]
        estimates = [_atom_marginal_correlator(probs_by_row[s.row]) for s in candidates]
        S[i, 0], S_err[i, 0] = _pool(estimates)

    # --- Photon marginals S_0j (j in {1,2,3}): pool rows with matching photon basis ---
    for photon, j in _PHOTON_BASIS_TO_INDEX.items():
        candidates = [s for s in rows_for_input if s.photon == photon]
        estimates = [
            _photon_marginal_correlator(probs_by_row[s.row], photon) for s in candidates
        ]
        S[0, j], S_err[0, j] = _pool(estimates)

    return S, S_err


# ---------------------------------------------------------------------------
# Density matrix assembly
# ---------------------------------------------------------------------------


def density_matrix_from_stokes(S: np.ndarray) -> np.ndarray:
    """rho = (1/4) sum_{ij} S_ij  sigma_i^A (x) sigma_j^P."""
    rho = np.zeros((4, 4), dtype=complex)
    for i in range(4):
        for j in range(4):
            rho += S[i, j] * np.kron(_PAULI[i], _PAULI[j])
    return rho / 4.0


def project_to_physical(rho: np.ndarray) -> np.ndarray:
    """
    Find the nearest physical (PSD, trace-1, Hermitian) state.

    Algorithm from Smolin, Gambetta, Smith (2012), PRL 108, 070502:
    eigen-decompose, clip negative eigenvalues to 0 with a uniform shift
    that preserves trace = 1.
    """
    rho = (rho + rho.conj().T) / 2  # enforce Hermiticity
    evals, evecs = np.linalg.eigh(rho)
    # Sort descending
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]

    n = len(evals)
    lam = evals.copy().real
    accumulator = 0.0
    i = n - 1
    while i >= 0 and lam[i] + accumulator / (i + 1) < 0:
        accumulator += lam[i]
        lam[i] = 0.0
        i -= 1
    for k in range(i + 1):
        lam[k] += accumulator / (i + 1)

    rho_phys = evecs @ np.diag(lam.astype(complex)) @ evecs.conj().T
    return (rho_phys + rho_phys.conj().T) / 2


# ---------------------------------------------------------------------------
# Target Bell states and fidelity
# ---------------------------------------------------------------------------


def _ket(coeffs: dict[str, complex]) -> np.ndarray:
    """Build a 4-dim ket in the |atom, photon> basis ordered as
    |0,V>, |0,pi>, |1,V>, |1,pi> — matching sigma_z eigenvalue convention
    (+1 first). coeffs keys: '0V', '0pi', '1V', '1pi'.
    """
    order = ["0V", "0pi", "1V", "1pi"]
    v = np.zeros(4, dtype=complex)
    for k, c in coeffs.items():
        v[order.index(k)] = c
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def _photon_R() -> tuple[complex, complex]:
    """|R> = (|V> + i|pi>)/sqrt(2)  in (V, pi) convention."""
    return (1 / np.sqrt(2), 1j / np.sqrt(2))


def _photon_L() -> tuple[complex, complex]:
    """|L> = (|V> - i|pi>)/sqrt(2)."""
    return (1 / np.sqrt(2), -1j / np.sqrt(2))


def target_bell(inp: Input) -> np.ndarray:
    """
    Build the ideal Bell state produced by CPF from |+>_atom (x) |R or L>_photon.

    |psi_R> = (1/sqrt(2)) (|1> |R>  -  |0> |L>)
    |psi_L> = (1/sqrt(2)) (|1> |L>  -  |0> |R>)
    """
    if inp == "R":
        r_V, r_pi = _photon_R()
        l_V, l_pi = _photon_L()
        # |1>|R> - |0>|L>
        coeffs = {
            "1V": r_V,
            "1pi": r_pi,
            "0V": -l_V,
            "0pi": -l_pi,
        }
    else:
        r_V, r_pi = _photon_R()
        l_V, l_pi = _photon_L()
        # |1>|L> - |0>|R>
        coeffs = {
            "1V": l_V,
            "1pi": l_pi,
            "0V": -r_V,
            "0pi": -r_pi,
        }
    # Normalize after the superposition is built
    # (the component amplitudes are already 1/sqrt(2) each)
    order = ["0V", "0pi", "1V", "1pi"]
    v = np.array([coeffs.get(k, 0) for k in order], dtype=complex)
    return v / np.linalg.norm(v)


def fidelity(rho: np.ndarray, psi: np.ndarray) -> float:
    """F = <psi|rho|psi> for a pure target psi."""
    return float(np.real(psi.conj() @ rho @ psi))


# ---------------------------------------------------------------------------
# Monte Carlo error propagation
# ---------------------------------------------------------------------------


def monte_carlo_fidelity(
    probs_by_row: dict[int, Probs],
    inp: Input,
    n_samples: int = 500,
    project_physical: bool = True,
    rng: np.random.Generator | None = None,
) -> tuple[float, float]:
    """
    Bootstrap fidelity by resampling each probability from a Gaussian with
    its reported uncertainty, renormalizing each row to 1, and recomputing.

    Returns (mean F, std F).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    psi = target_bell(inp)

    fids = []
    for _ in range(n_samples):
        sampled: dict[int, Probs] = {}
        for row, pr in probs_by_row.items():
            new_p: dict[tuple[str, str], tuple[float, float]] = {}
            # sample 4 joint probs independently, clip, then renormalize
            raw = {}
            for key, (p, e) in pr.p.items():
                s = max(0.0, rng.normal(p, e))
                raw[key] = s
            total = sum(raw.values())
            if total <= 0:
                # fall back to the central values
                new_p = pr.p.copy()
            else:
                for key, (p, e) in pr.p.items():
                    new_p[key] = (raw[key] / total, e)  # err not used further
            sampled[row] = Probs(new_p)
        S, _ = compute_stokes(sampled, inp)
        rho = density_matrix_from_stokes(S)
        if project_physical:
            rho = project_to_physical(rho)
        fids.append(fidelity(rho, psi))
    return float(np.mean(fids)), float(np.std(fids))


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------


def run_tomography(
    probs_by_row: dict[int, Probs],
    *,
    project_physical: bool = True,
    mc_samples: int = 500,
    verbose: bool = True,
) -> dict:
    """
    Full analysis: S matrices, density matrices, fidelities (with MC error)
    for both R and L inputs.
    """
    results: dict = {}

    for inp in ("R", "L"):
        S, S_err = compute_stokes(probs_by_row, inp)
        rho_linear = density_matrix_from_stokes(S)
        rho = project_to_physical(rho_linear) if project_physical else rho_linear
        psi = target_bell(inp)
        F_central = fidelity(rho, psi)
        F_mean, F_std = monte_carlo_fidelity(
            probs_by_row, inp, n_samples=mc_samples, project_physical=project_physical
        )

        results[inp] = {
            "S": S,
            "S_err": S_err,
            "rho_linear": rho_linear,
            "rho": rho,
            "target": psi,
            "F_central": F_central,
            "F_mc_mean": F_mean,
            "F_mc_std": F_std,
        }

        if verbose:
            print(f"\n=== Input: {inp} ===")
            print(f"S-matrix (rows = atom I/X/Y/Z, cols = photon I/X/Y/Z):")
            print(np.array2string(S, precision=3, suppress_small=True))
            print(f"Density matrix (after projection):")
            print(np.array2string(rho, precision=3, suppress_small=True))
            eigs = np.linalg.eigvalsh(rho).real
            print(f"Eigenvalues: {np.round(eigs, 4)}")
            print(f"Trace: {np.trace(rho).real:.4f}")
            print(f"Purity tr(rho^2): {np.trace(rho @ rho).real:.4f}")
            print(f"Fidelity (central): {F_central:.4f}")
            print(f"Fidelity (MC mean +/- std): {F_mean:.4f} +/- {F_std:.4f}")

    return results


# ---------------------------------------------------------------------------
# Integration helper: call your analyzer across all 18 rows
# ---------------------------------------------------------------------------


def collect_all_probs(
    analyzer,
    file_stem: str,
    params: dict,
    path: str = "./",
    post_select_sd: bool = False,
) -> dict[int, Probs]:
    """
    Call analyzer.reflection_analysis for each of the 18 rows and pack the
    results into a dict keyed by row number.

    file_stem example: '10_04_26_Bell_State_Prep_hail_marry_'
    The function appends '{row}' (no pickle suffix — your analyzer handles that).
    """
    probs_by_row: dict[int, Probs] = {}
    for row in range(1, 19):
        file_name = f"{file_stem}{row}"
        (
            P_4_atom_0,
            P_4_atom_0_err,
            P_4_atom_1,
            P_4_atom_1_err,
            P_7_atom_0,
            P_7_atom_0_err,
            P_7_atom_1,
            P_7_atom_1_err,
            _eff,
            _eff_err,
        ) = analyzer.reflection_analysis(
            file_name=file_name,
            path=path,
            parameters=params,
            plot_histogram=False,
            post_select_sd=post_select_sd,
        )
        probs_by_row[row] = Probs.from_reflection(
            P_4_atom_0,
            P_4_atom_0_err,
            P_4_atom_1,
            P_4_atom_1_err,
            P_7_atom_0,
            P_7_atom_0_err,
            P_7_atom_1,
            P_7_atom_1_err,
        )
    return probs_by_row


if __name__ == "__main__":
    analyzer = Analyzer(data_dir="./")

    ParamDictReflection_trap_off: ReflectionGateT = {
        "trigger_delay": 3.15e-6,
        "cooling_duration": 400e-6,
        "optical_pumping_duration": 340e-6,
        "pulse_delay": 19.2e-6,  # Trap on
        "pulse_duration": 1.2e-6,  # 1mus pulse
        "pulse_delay_SD_after": 21.5e-6,  # Trap on
        "pulse_duration_SD_after": 8e-6,
        "sequence_duration": 1.1e-3,
    }
    probs = collect_all_probs(
        analyzer,
        file_stem="10_04_26_Bell_State_Prep_hail_marry_",
        params=ParamDictReflection_trap_off,
    )
    results = run_tomography(probs)
    mle_results = run_mle(probs, mc_samples=100)

    # Access:
    rho_L_mle = mle_results["L"].rho
    F_L_mle = mle_results["L"].fidelity_mc_mean
    F_L_err = mle_results["L"].fidelity_mc_std

    inp = "L"
    rho = results[inp]["rho"]
    psi = results[inp]["target"]
    target = np.outer(psi, psi.conj())
    fig = plot_bell_state_tomography(
        rho_L_mle, target, inp, save_path=f"bell_state_{inp}_mle.pdf"
    )

    # for inp in ("R", "L"):
    #     rho = results[inp]["rho"]
    #     psi = results[inp]["target"]
    #     target = np.outer(psi, psi.conj())
    #     fig = plot_bell_state_tomography(
    #         rho, target, inp, save_path=f"bell_state_{inp}.svg"
    #     )
    #
    # plt.show()
